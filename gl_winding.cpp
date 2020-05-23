#include <memory>

#include "lbcommon.h"
#include <glm/gtc/type_ptr.hpp>


class WindingGLRenderer : public Renderer
{
public:
  // note that GLES does not support implicit type conversion, so, e.g. 1/2.0, x > 0, etc. doesn't work!
  //static constexpr
  const char* calcVS =
  R"(
  attribute vec2 normal_in;
  attribute vec2 va_in;
  attribute vec2 vb_in;

  uniform mat4 model;
  uniform float scale; // units per pixel
  uniform vec2 view_wh;
  uniform float y_min;

  varying vec2 va;
  varying vec2 vb;

  void main()
  {
    // screen coords for unexpanded vertices
    vec2 va_tf = vec2(model * vec4(va_in, 0.0f, 1.0f));
    vec2 vb_tf = vec2(model * vec4(vb_in, 0.0f, 1.0f));
    va = 0.5f*view_wh*(vec2(1.0f) + va_tf);
    vb = 0.5f*view_wh*(vec2(1.0f) + vb_tf);

    // normal.x > 0: we are right, otherwise left
    vec2 pos = (normal_in.x > 0.0f) == (va_tf.x > vb_tf.x) ? va_tf : vb_tf;

    // trapezoid will have less overdraw in some (most?) cases, but more complicated calc, so make a rectangle
    float y_max = max(va_tf.y, vb_tf.y);
    pos.y = normal_in.y > 0.0f ? y_max : y_min;

    // expand to get final position
    gl_Position = vec4(pos + 0.5f*(1.0f/400.0f)*normal_in, 0.0f, 1.0f);
  }
  )";

  const char* triangleCalcVS =
  R"(
  attribute vec2 position_in;
  attribute vec2 normal_in;
  attribute vec2 va_in;
  attribute vec2 vb_in;
  attribute vec2 vc_in;

  uniform mat4 model;
  uniform float scale; // units per pixel
  uniform vec2 view_wh;

  varying vec2 va;
  varying vec2 vb;
  varying vec2 vc;

  void main()
  {
    va = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(va_in, 0.0f, 1.0f)));
    vb = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vb_in, 0.0f, 1.0f)));
    vc = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vc_in, 0.0f, 1.0f)));

    // expand triangle by half a pixel so all partially covered pixels are included
    gl_Position = model * vec4(position_in + 0.5f*scale*normal_in, 0.0f, 1.0f);
  }
  )";

  const char* calcFS =
  R"(
  varying vec2 va;
  varying vec2 vb;
  //varying vec2 vc;

  // we expect v1.x >= v0.x and v0, v1 translated so pixel has corners (0,0) and (1,1)
  float areaEdge(vec2 v0, vec2 v1)
  {
    // edge entirely to left or right of pixel?
    if((v0.x < 0.0f && v1.x < 0.0f) || (v0.x > 1.0f && v1.x > 1.0f) || v0.x == v1.x)
      return 0.0f;

    vec2 dv = v1 - v0;
    float slope = dv.y/dv.x;
    // clip edge to pixel (x = 0 to 1)
    if(v0.x < 0.0f)
      v0 = vec2(0.0f, v0.y - slope*v0.x);
    if(v1.x > 1.0f)
      v1 = vec2(1.0f, v1.y + slope*(1.0f - v1.x));
    if(v0.y <= 0.0f && v1.y <= 0.0f)
      return 0.0f;
    if(v0.y >= 1.0f && v1.y >= 1.0f)
      return (v1.x - v0.x);  // *1, the height of the pixel
    // clip edge to bottom of pixel (y = 0)
    float invslope = dv.x/dv.y;  // 1/slope might be faster in GLSL
    if(v0.y < 0.0f)
      v0 = vec2(v0.x - invslope*v0.y, 0.0f);
    if(v1.y < 0.0f)
      v1 = vec2(v1.x - invslope*v1.y, 0.0f);
    if(v1.y > 1.0f) {
      float xi = v1.x + invslope*(1.0f - v1.y);
      return (v1.x - xi) + (xi - v0.x)*(v0.y + 1.0f)*0.5f;  // (v1.x - xi)*1 (height)
    }
    if(v0.y > 1.0f) {
      float xi = v0.x + invslope*(1.0f - v0.y);
      return (xi - v0.x) + (v1.x - xi)*(v1.y + 1.0f)*0.5f;  // (xi - v0.x)*1 (height)
    }
    // final case: clipped edge entirely inside triangle
    return (v1.x - v0.x)*(v1.y + v0.y)*0.5f;
  }

  /*
  float coverageTri(vec2 vA, vec2 vB, vec2 vC, vec2 p)
  {
    // translate so pixel has corners (0,0) and (1,1)
    p -= 0.5f;
    vA -= p;  vB -= p;  vC -= p;
    return (vA.x < vB.x ? areaEdge(vA, vB) : -areaEdge(vB, vA))
        + (vB.x < vC.x ? areaEdge(vB, vC) : -areaEdge(vC, vB))
        + (vC.x < vA.x ? areaEdge(vC, vA) : -areaEdge(vA, vC));
  }
  */

  void main()
  {
    // triangle winding
    //vec2 p = vec2(gl_FragCoord);
    //float coverage = coverageTri(va, vb, vc, p);
    //gl_FragColor = vec4(0.0f, 0.0f, 0.0f, coverage);

    // trapezoid winding
    vec2 p = vec2(gl_FragCoord);
    p -= 0.5f;
    float coverage = (va.x < vb.x ? areaEdge(va - p, vb - p) : -areaEdge(vb - p, va - p));
    gl_FragColor = vec4(coverage);
  }
  )";

  // fill (cover) pass
  const char* fillVS =
  R"(
  attribute vec2 position_in;

  uniform mat4 model;
  uniform float scale;
  uniform vec2 bbox_center;
  uniform vec2 bbox_wh;

  void main()
  {
    gl_Position = model * vec4(position_in*(0.5f*bbox_wh + vec2(scale)) + bbox_center, 0.0f, 1.0f);
  }
  )";

  const char* fillFS =
  R"(
  //#extension GL_EXT_shader_framebuffer_fetch : require

  uniform vec4 color;
  uniform vec2 view_wh;
  uniform sampler2D windingTex;
  uniform int mode;

  void main()
  {
    // all the actual work is done via blending
    if(mode == 0)
      gl_FragColor = color;
    else {
      // framebuffer texture (use .a instead of .r w/ glTextureBarrier())
      vec2 texcoord = gl_FragCoord.xy/view_wh;
      float W = texture2D(windingTex, texcoord).r;
      gl_FragColor = vec4(vec3(color), min(abs(W), 1.0f)*color.a);
    }

    // gl_LastFragData
    //float W = gl_LastFragData[0].a;
    //gl_FragColor = vec4(vec3(color), min(abs(W), 1.0f)*color.a);

    // debug
    //gl_FragColor = mix(vec4(1,0,0,1), vec4(0,0,1,1), 0.5*(W+1));
  }
  )";

  const char* screenVS =
  R"(
  attribute vec2 position_in;
  attribute vec2 texcoord_in;

  varying vec2 texcoord;

  void main()
  {
    texcoord = texcoord_in;
    gl_Position = vec4(position_in, 0.0, 1.0);
  }
  )";

  const char* screenFS =
  R"(
  uniform sampler2D texFramebuffer;

  varying vec2 texcoord;

  void main()
  {
    gl_FragColor = vec4(texture2D(texFramebuffer, texcoord).rgb, 1.0);

    // debug
    //float W = texture2D(texFramebuffer, texcoord).a;
    //gl_FragColor = mix(vec4(1,0,0,1), vec4(0,0,1,1), 0.5*(W+1));

    ////vec2 texcoord2 = gl_FragCoord.xy/vec2(800, 800);
    //float W = texture2D(texFramebuffer, texcoord).a;
    //gl_FragColor = vec4(0, 1, 0, min(abs(W), 1.0f));
  }
  )";

  std::unique_ptr<ShaderProgram> spCalc;
  std::unique_ptr<ShaderProgram> spFill;
  std::unique_ptr<ShaderProgram> spScreen;

  GLuint frameBuffer;
  GLuint texFBWinding;
  struct { GLuint position, normal, va, vb, vc; } attrCalc;
  struct { GLint model, scale, view_wh, y_min; } unifCalc;
  struct { GLuint position; } attrFill;
  struct { GLint model, scale, view_wh, bbox_center, bbox_wh, mode, color; } unifFill;

  WindingGLRenderer()
  {
    name = "GL Winding";
    // shader setup
    spCalc.reset(new ShaderProgram(calcVS, calcFS));
    attrCalc = { spCalc->attrib("position_in"), spCalc->attrib("normal_in"),
        spCalc->attrib("va_in"), spCalc->attrib("vb_in"), spCalc->attrib("vc_in") };
    unifCalc = { spCalc->uniform("model"),
        spCalc->uniform("scale"), spCalc->uniform("view_wh"), spCalc->uniform("y_min") };
    // glEnableVertexAttribArray actually applies to the VAO, not shader, so calling for each shader attrib
    //  w/o a (different) VAO doesn't make sense ... but do it anyway for now
    glEnableVertexAttribArray(attrCalc.normal);
    glEnableVertexAttribArray(attrCalc.va);
    glEnableVertexAttribArray(attrCalc.vb);
    // for triangle winding renderer
    //glEnableVertexAttribArray(attrCalc.position);
    //glEnableVertexAttribArray(attrCalc.vc);

    spFill.reset(new ShaderProgram(fillVS, fillFS));
    attrFill = { spFill->attrib("position_in") };
    unifFill = { spFill->uniform("model"), spFill->uniform("scale"), spFill->uniform("view_wh"),
        spFill->uniform("bbox_center"), spFill->uniform("bbox_wh"), spFill->uniform("mode"), spFill->uniform("color") };
    glEnableVertexAttribArray(attrFill.position);
    glUniform1i(spFill->uniform("windingTex"), 0);

    spScreen.reset(new ShaderProgram(screenVS, screenFS));
    glEnableVertexAttribArray(spScreen->attrib("position_in"));
    glEnableVertexAttribArray(spScreen->attrib("texcoord_in"));

    // framebuffer setup
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    // Create texture for winding buffer
    glGenTextures(1, &texFBWinding);
    glBindTexture(GL_TEXTURE_2D, texFBWinding);
    // offically, blending only works for 16F, not 32F, but 32F seems to work (although >2x slower)
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, 800, 800, 0, GL_RGBA, GL_HALF_FLOAT, NULL);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, 800, 800, 0, GL_RED, GL_HALF_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texFBWinding, 0);

    // Create Renderbuffer Object to hold depth and stencil buffers
    //GLuint rboDepthStencil;
    //glGenRenderbuffers(1, &rboDepthStencil);
    //glBindRenderbuffer(GL_RENDERBUFFER, rboDepthStencil);
    //glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 800, 800);
    //glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rboDepthStencil);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
    glDisable(GL_DEPTH_TEST);
    //glDepthMask(GL_FALSE);  -- disabling depth test apparently also disables writes

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      PLATFORM_LOG("Error creating framebuffer\n");
    GL_CHECK(glBindFramebuffer(GL_FRAMEBUFFER, 0));
  }

  ~WindingGLRenderer() override
  {
    glDeleteTextures(1, &texFBWinding);
    //glDeleteTextures(1, &texFBColor);
    //glDeleteRenderbuffers
    glDeleteFramebuffers(1, &frameBuffer);
  }

  void beginFrame() override
  {
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
    //glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT); // | GL_DEPTH_BUFFER_BIT);
  }

  void drawGlyph(LBGlyph* g) override
  {
    // draw triangles
    glUseProgram(spCalc->prog);
    glBlendEquation(GL_FUNC_ADD);
    //glBlendFuncSeparate(GL_ZERO, GL_ONE, GL_ONE, GL_ONE);  // dest.rgb = 0*src.rgb + 1*dest.rgb
    glBlendFunc(GL_ONE, GL_ONE);

    // note that model matrix is different for each glyph
    glUniformMatrix4fv(unifCalc.model, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform1f(unifCalc.scale, UniformBase::scale);
    glUniform2f(unifCalc.view_wh, UniformBase::view_wh.x, UniformBase::view_wh.y);

    // trapezoid (rectangle) winding calc
    // bottom of trapezoids in transformed coords
    vec4 c0 = UniformBase::model * vec4(g->x_min, g->y_min, 0, 1);
    vec4 c1 = UniformBase::model * vec4(g->x_max, g->y_min, 0, 1);
    vec4 c2 = UniformBase::model * vec4(g->x_min, g->y_max, 0, 1);
    vec4 c3 = UniformBase::model * vec4(g->x_max, g->y_max, 0, 1);
    float y_min = glm::min(glm::min(c0.y, c1.y), glm::min(c2.y, c3.y));
    glUniform1f(unifCalc.y_min, y_min);

    // vertices
    glVertexAttribPointer(attrCalc.normal, 2, GL_FLOAT, GL_FALSE, sizeof(TrapezoidVertex), &g->trapezoid_vertices[0].triangle_normal);
    glVertexAttribPointer(attrCalc.va, 2, GL_FLOAT, GL_FALSE, sizeof(TrapezoidVertex), &g->trapezoid_vertices[0].va);
    glVertexAttribPointer(attrCalc.vb, 2, GL_FLOAT, GL_FALSE, sizeof(TrapezoidVertex), &g->trapezoid_vertices[0].vb);
    glDrawArrays(GL_TRIANGLES, 0, g->num_trapezoid_vertices);

    /*
    // triangle winding calc
    glVertexAttribPointer(attrCalc.position, 2, LB_REAL, GL_FALSE, sizeof(Vertex), &g->winding_vertices[0].x);
    glVertexAttribPointer(attrCalc.normal, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->winding_vertices[0].edge_normal);
    glVertexAttribPointer(attrCalc.va, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->winding_vertices[0].va);
    glVertexAttribPointer(attrCalc.vb, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->winding_vertices[0].vb);
    glVertexAttribPointer(attrCalc.vc, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->winding_vertices[0].vc);
    glDrawArrays(GL_TRIANGLES, 0, g->num_winding_vertices);
    */
    /*
    // tessellated triangles
    glVertexAttribPointer(attrCalc.position, 2, LB_REAL, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].x);
    glVertexAttribPointer(attrCalc.normal, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].edge_normal);
    glVertexAttribPointer(attrCalc.va, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].va);
    glVertexAttribPointer(attrCalc.vb, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].vb);
    glVertexAttribPointer(attrCalc.vc, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), &g->triangle_vertices[0].vc);
    glDrawArrays(GL_TRIANGLES, 0, g->num_triangle_vertices);
    */

    // original idea was to use two passes of blending to compute min(abs(W), 1)), w/ abs(W) = max(W, -W),
    //  but for some absurd reason, GL_MIN and GL_MAX blend equations do not use glBlendFunc factors
    // Fortunately, it looks like swapping framebuffer gives basically same performance as gl_LastFragData

    glUseProgram(spFill->prog);
    glUniformMatrix4fv(unifFill.model, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform1f(unifFill.scale, UniformBase::scale);
    glUniform2f(unifFill.view_wh, UniformBase::view_wh.x, UniformBase::view_wh.y);
    glUniform2f(unifFill.bbox_center, (g->x_max + g->x_min)/2, (g->y_max + g->y_min)/2);
    glUniform2f(unifFill.bbox_wh, g->x_max - g->x_min, g->y_max - g->y_min);

    glVertexAttribPointer(attrFill.position, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][0]);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texFBWinding);
    glUniform1i(unifFill.mode, 1); // winding to coverage

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glUniform4fv(unifFill.color, 1, glm::value_ptr(UniformBase::color));
    glDrawArrays(GL_TRIANGLES, 0, 6);

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    // reset winding buffer ... performance seems OK (were it not, we could try glScissor + glClear)
    glUniform1i(unifFill.mode, 0); // color pass-thru (clear)
    glUniform4f(unifFill.color, 0, 0, 0, 1);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    //glClear(GL_COLOR_BUFFER_BIT);  // slower

    /*
    // single pass fill for use w/ glTextureBarrier, or better, gl_LastFragData
    // if color.a = 1 everywhere and shape is simple (single contour, no self-intersections), we can use
    //  `GL_DST_ALPHA ...` to apply coverage from alpha channel and clear in a single pass
    //GL_CHECK(glTextureBarrier());
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ZERO, GL_ZERO);
    //glBlendFuncSeparate(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, GL_ZERO, GL_ZERO);
    glUniform4fv(unifFill.color, 1, glm::value_ptr(UniformBase::color));
    glDrawArrays(GL_TRIANGLES, 0, 6);
    */
  }

  void endFrame() override
  {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);  // draw to screen
    /*
    // reset blend state; we could disable GL_BLEND for this pass, but for now we use gl_FragColor.a = 1 in FS
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(spScreen->prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texFBWinding);
    glUniform1i(spScreen->uniform("texFramebuffer"), 0);

    glVertexAttribPointer(spScreen->attrib("position_in"), 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][0]);
    glVertexAttribPointer(spScreen->attrib("texcoord_in"), 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][2]);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    // unnecessary?
    glBindTexture(GL_TEXTURE_2D, 0);
    */
  }
};
