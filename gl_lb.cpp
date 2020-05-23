#include <memory>

#include "lbcommon.h"
#include "glpp.h"
#include <glm/gtc/type_ptr.hpp>


class LBGLRenderer : public Renderer
{
public:
  //static constexpr
  const char* fillVS =
  R"(#version 120

  attribute vec2 position;
  attribute vec2 normal;
  attribute vec3 fringe;

  uniform mat4 model;
  uniform float scale; // units per pixel

  varying vec3 edge;
  varying vec3 edgescale;

  void main()
  {
    float s2 = model[2][2];
    edge = fringe*s2;
    edgescale = abs(fringe*s2); // + vec3(fringe == vec3(0.0));
    // expand triangle by half a pixel so all partially covered pixels are included
    gl_Position = model * vec4(position + 0.5*scale*normal, 0.0, 1.0);
  }
  )";

  const char* fillFS =
  R"(#version 120

  uniform vec4 color;
  uniform float aastep;
  varying vec3 edge;
  varying vec3 edgescale;

  void main()
  {
    // remember that value of edge here is interpolated from 3 vertices defining our triangle
    // note that we could shift this calculation to vs, letting interpolation generate alpha directly
    vec3 vv = edgescale - edge; //abs(edge);
//    float v = min(vv.x, min(vv.y, vv.z));
//    float a = v > 0 ? smoothstep(0.0f, aastep, v) : 1.0f;

    float v = min(vv.x > 0 ? vv.x : 1.0f, min(vv.y > 0 ? vv.y : 1.0f, vv.z > 0 ? vv.z : 1.0f));
    float a = v > 0 ? smoothstep(0.0f, aastep, v) : 1.0f;

    gl_FragColor = vec4(color.rgb, a*color.a);
    // this is handy for debugging
    //gl_FragColor = mix(vec4(2.0, 0.0, 0.0, 1.0), vec4(0.0, 0.0, 2.0, 1.0), clamp(0.5 + 0.5*edge.x/edgescale.x, 0, 1));
  }
  )";

  const char* lbVS =
  R"(#version 120

  attribute vec3 position;
  attribute vec3 texcoord;
  uniform mat4 model;
  varying vec3 texCoord;

  void main()
  {
    texCoord = texcoord;
    gl_Position = model * vec4(position.xy, 0.0, 1.0);
  }
  )";

  const char* lbFS =
  R"(#version 120

  uniform vec4 color;
  varying vec3 texCoord;

  void main()
  {
    vec2 p, px, py;
    float fx,fy,sd,alpha;

    p = texCoord.xy;
    px = dFdx(p);
    py = dFdy(p);
    fx = 2.0*p.x*px.x - px.y;
    fy = 2.0*p.x*py.x - py.y;
    sd = (p.x*p.x  - p.y) / sqrt(fx*fx + fy*fy);
    alpha = 0.5 - sd * texCoord.z;
    alpha = clamp(alpha, 0.0, 1.0);

    //frag_color = vec4(p.x,p.y,0.0,1.0);
    gl_FragColor = vec4(color.rgb, alpha*color.a);
  }
  )";

  std::unique_ptr<ShaderProgram> spFill;
  GLuint attrPos_Fill, attrFringe_Fill, attrNormal_Fill;
  GLint uniModel_Fill, uniScale_Fill, uniAAStep_Fill, uniColor_Fill;

  std::unique_ptr<ShaderProgram> spLB;
  GLuint attrPos_LB, attrFringe_LB;
  GLint uniModel_LB, uniColor_FB;

  // TODO: maybe have ShaderProgram cache the uniform, attrib indices (unordered_map?) so we don't have to
  LBGLRenderer()
  {
    spFill.reset(new ShaderProgram(fillVS, fillFS));
    attrPos_Fill = spFill->attrib("position");
    glEnableVertexAttribArray(attrPos_Fill);
    // If an attribute or uniform is not used in shader, it may be optimized away and trying to access it will fail
    attrFringe_Fill = spFill->attrib("fringe");
    glEnableVertexAttribArray(attrFringe_Fill);
    attrNormal_Fill = spFill->attrib("normal");
    glEnableVertexAttribArray(attrNormal_Fill);
    // uniforms - must call glUseProgram before any uniforms can be set!
    uniModel_Fill = spFill->uniform("model");  // VS
    uniScale_Fill = spFill->uniform("scale");  // VS
    uniAAStep_Fill = spFill->uniform("aastep");  // FS
    uniColor_Fill = spFill->uniform("color");

    spLB.reset(new ShaderProgram(lbVS, lbFS));
    attrPos_LB = spLB->attrib("position");
    glEnableVertexAttribArray(attrPos_LB);
    // If an attribute or uniform is not used in shader, it may be optimized away and trying to access it will fail
    attrFringe_LB = spLB->attrib("texcoord");
    glEnableVertexAttribArray(attrFringe_LB);
    // uniforms
    uniModel_LB = spLB->uniform("model");
    uniColor_FB = spLB->uniform("color");
  }

  void drawGlyph(LBGlyph* g) override
  {
    // draw triangles
    spFill->use();
    // uniforms
    glUniformMatrix4fv(uniModel_Fill, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform1f(uniScale_Fill, UniformBase::scale);
    glUniform1f(uniAAStep_Fill, UniformBase::aastep);
    glUniform4fv(uniColor_Fill, 1, glm::value_ptr(UniformBase::color));
    // vertices
    glVertexAttribPointer(attrPos_Fill, 2, LB_REAL, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].x);
    glVertexAttribPointer(attrFringe_Fill, 3, LB_REAL, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].u);
    glVertexAttribPointer(attrNormal_Fill, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].triangle_normal);
    glDrawArrays(GL_TRIANGLES, 0, g->num_triangle_vertices);

    // draw inside curves
    spLB->use();
    // uniforms
    glUniformMatrix4fv(uniModel_LB, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform4fv(uniColor_FB, 1, glm::value_ptr(UniformBase::color));
    // vertices
    glVertexAttribPointer(attrPos_LB, 3, LB_REAL, GL_FALSE, sizeof(Vertex), &g->inside_curve_triangles[0].x);
    glVertexAttribPointer(attrFringe_LB, 3, LB_REAL, GL_FALSE, sizeof(Vertex), &g->inside_curve_triangles[0].u);
    glDrawArrays(GL_TRIANGLES, 0, g->num_inside_ctrl_points*3);

    // draw outside curves
    glVertexAttribPointer(attrPos_LB, 3, LB_REAL, GL_FALSE, sizeof(Vertex), &g->outside_curve_triangles[0].x);
    glVertexAttribPointer(attrFringe_LB, 3, LB_REAL, GL_FALSE, sizeof(Vertex), &g->outside_curve_triangles[0].u);
    glDrawArrays(GL_TRIANGLES, 0, g->num_outside_ctrl_points*3);

    // draw some stuff for debugging
#if 0
    // draw triangle outlines
    /* glColor3f(1,1,1);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, g->num_triangle_indices , GL_UNSIGNED_INT, g->triangle_indices);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    */
    // draw the control points
    drawGlyphData(g);
#endif
  }
};

/*
class SDFGLRenderer : Renderer
{
public:
  // identical to fillVS
  const char* fillSDFVS =
  R"(#version 120

  attribute vec2 position;
  attribute vec2 normal;
  attribute vec3 fringe;

  uniform mat4 model;
  uniform float scale; // units per pixel

  varying vec3 edge;
  varying vec3 edgescale;

  void main()
  {
    float s2 = model[2][2];
    edge = fringe*s2;
    edgescale = abs(fringe*s2); // + vec3(fringe == vec3(0.0));
    // expand triangle by half a pixel so all partially covered pixels are included
    gl_Position = model * vec4(position + 0.5*scale*normal, 0.0, 1.0);
  }
  )";

  // same as fillFS, except it writes out the distance directly; to be used with glBlendEquation(GL_MAX)
  const char* fillSDFGenFS =
  R"(#version 120

  uniform vec4 color;
  uniform float aastep;
  varying vec3 edge;
  varying vec3 edgescale;

  void main()
  {
    // remember that value of edge here is interpolated from 3 vertices defining our triangle
    vec3 vv = edgescale - edge;
    float v = min(vv.x > 0 ? vv.x : 1.0f, min(vv.y > 0 ? vv.y : 1.0f, vv.z > 0 ? vv.z : 1.0f));
    gl_FragColor = vec4(v, 0, 0, 1.0f);
  }
  )";

  // either use same geometry for second pass (SDF), or have third pass to clear texture if using bbox quad
  const char* sdfVS =
  R"(#version 120

  attribute vec2 position;

  attribute vec2 normal;

  attribute vec2 texcoord;

  uniform mat4 model;
  uniform float scale; // units per pixel

  varying vec2 texCoord;

  void main()
  {
    texCoord = texcoord;
    //gl_Position = vec4(position, 0.0, 1.0);
    gl_Position = model * vec4(position + 0.5*scale*normal, 0.0, 1.0);
  }
  )";

  const char* sdfFS =
  R"(#version 120

  uniform vec4 color;
  uniform float aastep;
  uniform sampler2D texSDF;
  varying vec2 texCoord;

  void main()
  {
    float v = texture2D(texSDF, texCoord).r;
    float a = v > 0 ? smoothstep(0.0f, aastep, v) : 1.0f;
    gl_FragColor = vec4(color.rgb, a*color.a);
  }
  )";

};
*/
