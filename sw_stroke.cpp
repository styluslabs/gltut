#include "swrender.h"

namespace Stroke {

using namespace glm;

struct Varying : public VaryingBase {
  float edge;
  float edgescale;
};

// putting uniforms outside shader classes mirrors the fact that GL VS and FS uniforms names share the same namespace
const Vertex* verts;

// VS uniforms
float strokeWidth;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    // inputs
    position((float)verts[idx].x, (float)verts[idx].y), normal(verts[idx].edge_normal), direction((float)verts[idx].u),
    // outputs
    edge(vout->edge), edgescale(vout->edgescale) {}

  // in
  vec2 position;
  const vec2& normal;
  float direction;

  // out
  float& edge;
  float& edgescale;

  void main()
  {
    float s2 = model[2][2];
    float extrude = direction*(0.5f + strokeWidth*s2);
    edge = extrude;
    edgescale = abs(extrude);
    gl_Position = model * vec4(position + extrude*scale*normal, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    edge(baryinterp(vA.edge, vB.edge, vC.edge, t)),
    edgescale(baryinterp(vA.edgescale, vB.edgescale, vC.edgescale, t)) {}

  // in
  float edge;
  float edgescale;

  void main()
  {
    float vv = edgescale - abs(edge);
    float a = clamp(vv, 0.0f, 1.0f);  //smoothstep(0.0f, 1.0f, vv);
    gl_FragColor = vec4(vec3(color), a*color.a);
  }
};

} // end namespace Stroke


namespace Stroke_LB {

using namespace glm;

struct Varying : public VaryingBase {
  vec3 texCoord;
};

const Vertex* verts;

// VS uniforms
float strokeWidth;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position((float)verts[idx].x, (float)verts[idx].y, (float)verts[idx].z),
    texcoord((float)verts[idx].u, (float)verts[idx].v, (float)verts[idx].w),
    normal(verts[idx].triangle_normal),
    texCoord(vout->texCoord) {}

  // in
  vec3 position;
  vec3 texcoord;
  const vec2& normal;

  // out
  vec3& texCoord;

  void main()
  {
    // TODO: expand triangle, scale texCoord accordingly - switch(texCoord.x) 0: 0.5: 1:? or figure out more general formula?
    float s2 = model[2][2];
    float extrude = 0.5f + strokeWidth*s2;

    texCoord = texcoord;
    gl_Position = model * vec4(vec2(position) + texcoord.z*extrude*scale*normal, 0.0, 1.0);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)),
    dFdx_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.x, t.db.x, t.dc.x)),
    dFdy_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.y, t.db.y, t.dc.y)) {}

  // in
  vec3 texCoord;
  vec3 dFdx_texCoord;
  vec3 dFdy_texCoord;

  void main()
  {
    vec2 p = vec2(texCoord);
    if(p.x < 0 || p.x > 1 || p.y < 0 || p.y > 1)
      return;
    vec2 px = vec2(dFdx_texCoord);  // dFdx(p);
    vec2 py = vec2(dFdy_texCoord);  // dFdy(p);
    float fx = 2.0f*p.x*px.x - px.y;
    float fy = 2.0f*p.x*py.x - py.y;
    float sd = abs((p.x*p.x  - p.y) / sqrt(fx*fx + fy*fy)) - strokeWidth*model[2][2];
    float alpha = clamp(0.5f - sd, 0.0f, 1.0f);
    gl_FragColor = vec4(vec3(color), alpha*color.a);
  }
};

} // end namespace Stroke_LB


// render a circle
namespace RoundJoin {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 xy;
  //float edgescale;
};

// putting uniforms outside shader classes mirrors the fact that GL VS and FS uniforms names share the same namespace
const Vertex* verts;

// VS uniforms
float strokeWidth;  // actually HALF width

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    // inputs
    position((float)verts[idx].x, (float)verts[idx].y), direction((float)verts[idx].u, (float)verts[idx].v),
    // outputs
    xy(vout->xy) {} //, edgescale(vout->edgescale) {}

  // in
  vec2 position;
  vec2 direction;

  // out
  vec2& xy;
  //float& edgescale;

  void main()
  {
    float s2 = model[2][2];
    vec2 extrude = direction*(0.5f + strokeWidth*s2);
    xy = extrude;
    //edgescale = abs(extrude);
    gl_Position = model * vec4(position + scale*extrude, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    xy(baryinterp(vA.xy, vB.xy, vC.xy, t)) {}

  // in
  vec2 xy;
  //float edgescale;

  void main()
  {
    // in the future, we could also calculate the angle of the current point to limit rendering to an arc
    float d = length(xy) - strokeWidth*model[2][2];  // dist from edge of circle
    float a = 0.5f - clamp(d, -0.5f, 0.5f);
    gl_FragColor = vec4(vec3(color), a*color.a);
  }
};

} // end namespace RoundJoin


class StrokeSWRenderer : public SWRenderer
{
public:
  // drawElements indices to generate triangles for quads
  std::vector<int> strokeElems;

  StrokeSWRenderer()
  {
    //strokeElems.reserve(6000);
    for(int ii = 0; ii < 1000; ii += 4) {
      strokeElems.push_back(ii + 0);
      strokeElems.push_back(ii + 1);
      strokeElems.push_back(ii + 2);
      strokeElems.push_back(ii + 2);
      strokeElems.push_back(ii + 3);
      strokeElems.push_back(ii + 0);
    }
  }

  void drawGlyph(LBGlyph* g) override
  {
    // uniforms
    Stroke::strokeWidth = 6.0f;  // actually the half width
    // vertices
    Stroke::verts = g->stroke_vertices + 2;

    RoundJoin::strokeWidth = 6.0f;  // actually the half width
    RoundJoin::verts = g->stroke_vertices;
    ShaderPipeline<Stroke::VS, Stroke::FS, Stroke::Varying, BlendOver> StrokePipeline(blendOver);
    ShaderPipeline<RoundJoin::VS, RoundJoin::FS, RoundJoin::Varying, BlendOver> JoinPipeline(blendOver);
    for(int ii = 0; ii < g->num_contours; ++ii) {
      int nelems = ((ii+1 < g->num_contours) ? g->contours[ii+1] : g->num_verts) - g->contours[ii];
      // draw segments
      drawElements(StrokePipeline, &strokeElems[0], 6*nelems, 4*nelems);
      // draw joins - this actually draws starting (= ending) join twice
      drawElements(JoinPipeline, &strokeElems[0], 6*(nelems + 1), 4*(nelems + 1));
      // advance to next contour
      Stroke::verts += 4*(nelems + 1);
      RoundJoin::verts += 4*(nelems + 1);
    }

    // draw curved segments
    Stroke_LB::strokeWidth = 6.0f;  // actually the half width
    ShaderPipeline<Stroke_LB::VS, Stroke_LB::FS, Stroke_LB::Varying, BlendOver> LBPipeline(blendOver);
    // inside curves
    Stroke_LB::verts = g->inside_curve_triangles;
    drawArrays(LBPipeline, g->num_inside_ctrl_points*3);
    // outside curves
    Stroke_LB::verts = g->outside_curve_triangles;
    drawArrays(LBPipeline, g->num_outside_ctrl_points*3);
  }
};
