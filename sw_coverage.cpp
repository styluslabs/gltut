#include "swrender.h"
#include <vector>

class CoverageAccum {
public:
  Buffer2D<float>& coverageBuffer;
  CoverageAccum(Buffer2D<float>& buff) : coverageBuffer(buff) {}

  // we expect FS to write distance to gl_FragColor.a
  void blend(int x, int y, const FSBase& fs)
  {
    float& dest = coverageBuffer.data[x + coverageBuffer.width*y];
    dest = dest + fs.gl_FragColor.a;
  }
};

// Sutherland-Hodgman polygon clipping algorithm - clipping polygon must be convex
//  ref: http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Lua
// Keep in mind that we need to handle the case of arbitrarily small triangles - we can't assume triangle
//  covers several fragments

static bool halfPlane(const vec2& p, const vec2& cp1, const vec2& cp2)
{
  return (cp2.x-cp1.x)*(p.y-cp1.y) > (cp2.y-cp1.y)*(p.x-cp1.x);
}

static vec2 intersection(const vec2& cp1, const vec2& cp2, const vec2& s, const vec2& e)
{
  vec2 dc = cp1 - cp2;
  vec2 dp = s - e;
  float n1 = cp1.x*cp2.y - cp1.y*cp2.x;
  float n2 = s.x*e.y - s.y*e.x;
  float n3 = 1 / (dc.x*dp.y - dc.y*dp.x);
  float x = (n1*dp.x - n2*dc.x) * n3;
  float y = (n1*dp.y - n2*dc.y) * n3;
  return vec2(x,y);
}

// cp1.y == cp2.y
//static vec2 intersectHorz(const vec2& cp1, const vec2& cp2, const vec2& s, const vec2& e)
//{
//  float x = (cp1.y - e.y)*(s.x - e.x)/(s.y - e.y) + e.x;
//  return vec2(x, cp1.y);
//}

static std::vector<vec2> clipPolygon(const std::vector<vec2>& subject, const std::vector<vec2>& clip)
{
  std::vector<vec2> out = subject;
  // iterate over edges of clip polygon
  for(size_t cp1 = clip.size() - 1, cp2 = 0; cp2 < clip.size(); cp1 = cp2++) {
    std::vector<vec2> in = out;
    out.clear();
    // iterate over edges of subject polygon
    for(size_t e = 0, s = in.size() - 1; e < in.size(); s = e++) {
      if(halfPlane(in[e], clip[cp1], clip[cp2])) {
        if(!halfPlane(in[s], clip[cp1], clip[cp2]))
          out.push_back(intersection(clip[cp1], clip[cp2], in[s], in[e]));
        out.push_back(in[e]);
      }
      else if(halfPlane(in[s], clip[cp1], clip[cp2]))
        out.push_back(intersection(clip[cp1], clip[cp2], in[s], in[e]));
    }
  }
  return out;
}

// calculate signed area of any polygon - http://alienryderflex.com/polygon_area/
// also: https://en.wikipedia.org/wiki/Shoelace_formula (a special case of Green's theorem)
static float polygonArea(const std::vector<vec2>& points)
{
  float area = 0;
  for(size_t ii = 0, jj = points.size() - 1; ii < points.size(); jj = ii++)
    area += (points[jj].x + points[ii].x)*(points[jj].y - points[ii].y);
  return area/2;
}

// optimized calculation of triangle - pixel square intersection area; sums signed area between each edge
//  and x axis which falls inside pixel square
// some inspiration from http://www.cap-lore.com/MathPhys/IP/ (general polygon intersection area)
// other options might be "pipelined" version of Sutherland-Hodgeman w/ optimizations for convex-convex case,
//  or the linear-time convex-convex algorithm from Computational Geometry in C by O'Rourke
// for GLSL, we should try to eliminate as many branches as possible

// we expect v1.x >= v0.x and v0, v1 translated so pixel has corners (0,0) and (1,1)
float areaEdge(vec2 v0, vec2 v1)
{
  // edge entirely to left or right of pixel?
  if((v0.x < 0 && v1.x < 0) || (v0.x > 1 && v1.x > 1) || v0.x == v1.x)
    return 0.0f;

  vec2 dv = v1 - v0;
  float slope = dv.y/dv.x;
  // clip edge to pixel (x = 0 to 1)
  if(v0.x < 0)
    v0 = vec2(0.0f, v0.y - slope*v0.x);
  if(v1.x > 1)
    v1 = vec2(1.0f, v1.y + slope*(1.0f - v1.x));
  if(v0.y <= 0 && v1.y <= 0)
    return 0.0f;
  if(v0.y >= 1 && v1.y >= 1)
    return (v1.x - v0.x);  // *1, the height of the pixel
  // clip edge to bottom of pixel (y = 0)
  float invslope = dv.x/dv.y;  // 1/slope might be faster in GLSL
  if(v0.y < 0)
    v0 = vec2(v0.x - invslope*v0.y, 0.0f);
  if(v1.y < 0)
    v1 = vec2(v1.x - invslope*v1.y, 0.0f);
  if(v1.y > 1) {
    float xi = v1.x + invslope*(1.0f - v1.y);
    return (v1.x - xi) + (xi - v0.x)*(v0.y + 1.0f)/2;  // (v1.x - xi)*1 (height)
  }
  if(v0.y > 1) {
    float xi = v0.x + invslope*(1.0f - v0.y);
    return (xi - v0.x) + (v1.x - xi)*(v1.y + 1.0f)/2;  // (xi - v0.x)*1 (height)
  }
  // final case: clipped edge entirely inside triangle
  return (v1.x - v0.x)*(v1.y + v0.y)/2;
}

// not static - used in sw_winding.cpp
float coverageTri(vec2 vA, vec2 vB, vec2 vC, vec2 p)
{
  // translate so pixel has corners (0,0) and (1,1)
  p -= 0.5f;
  vA -= p;  vB -= p;  vC -= p;
  return (vA.x < vB.x ? areaEdge(vA, vB) : -areaEdge(vB, vA))
      + (vB.x < vC.x ? areaEdge(vB, vC) : -areaEdge(vC, vB))
      + (vC.x < vA.x ? areaEdge(vC, vA) : -areaEdge(vA, vC));
}

namespace CoverageCalc {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 va, vb, vc;
};

const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
  // edge_normal expands triangle in all directions
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position_in((float)verts[idx].x, (float)verts[idx].y), normal_in(verts[idx].edge_normal),
    va_in(verts[idx].va), vb_in(verts[idx].vb), vc_in(verts[idx].vc),
    va(vout->va), vb(vout->vb), vc(vout->vc) {}

  // in
  vec2 position_in;
  const vec2& normal_in;
  const vec2& va_in;
  const vec2& vb_in;
  const vec2& vc_in;

  // out
  vec2& va;
  vec2& vb;
  vec2& vc;

  void main()
  {
    // screen coords for unexpanded vertices
    va = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(va_in, 0.0f, 1.0f)));
    vb = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vb_in, 0.0f, 1.0f)));
    vc = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vc_in, 0.0f, 1.0f)));

    // expand triangle by half a pixel so all partially covered pixels are included
    gl_Position = model * vec4(position_in + 0.5f*scale*normal_in, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FrontFacing(cross(vec2(vB.position) - vec2(vA.position), vec2(vC.position) - vec2(vA.position)) > 0),
    gl_FragCoord(t.x + 0.5f, t.y + 0.5f, 0.0f, 1.0f),
    va(baryinterp(vA.va, vB.va, vC.va, t)),
    vb(baryinterp(vA.vb, vB.vb, vC.vb, t)),
    vc(baryinterp(vA.vc, vB.vc, vC.vc, t)) {}

  bool gl_FrontFacing;
  vec4 gl_FragCoord;

  // in
  vec2 va;
  vec2 vb;
  vec2 vc;

  void main()
  {
    vec2 p = vec2(gl_FragCoord);
    float coverage = coverageTri(va, vb, vc, p);
    gl_FragColor = vec4(0, 0, 0, coverage);

    //std::vector<vec2> clipped = clipPolygon({va, vb, vc},
    //    {vec2(p.x-0.5f,p.y-0.5f), vec2(p.x+0.5f,p.y-0.5f), vec2(p.x+0.5f,p.y+0.5f), vec2(p.x-0.5f,p.y+0.5f)});
    //float coverage0 = clipped.size() > 2 ? polygonArea(clipped) : 0.0f;
    //if(abs(coverage - coverage0) > 0.1f)
    //  fprintf(stderr, "coverage0 = %f, coverage = %f\n", coverage0, coverage);
  }
};

} // end namespace CoverageCalc

namespace CoverageFill {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 texCoord;
};

// vertex: {x,y,u,v} - this shader will typically just be used with a pair of triangles forming a quad
const float (*verts)[4];

// additional uniforms
vec2 bbox_center;
vec2 bbox_wh;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position_in(verts[idx][0], verts[idx][1]), texcoord_in(verts[idx][2], verts[idx][3]),
    texCoord(vout->texCoord) {}

  // in
  vec2 position_in;
  vec2 texcoord_in;

  // out
  vec2& texCoord;

  void main()
  {
    texCoord = texcoord_in;
    gl_Position = model * vec4(position_in*(0.5f*bbox_wh + vec2(scale)) + bbox_center, 0.0f, 1.0f);
  }
};

// FS uniforms
Buffer2D<float> coverageTex;

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FragCoord(t.x + 0.5, t.y + 0.5), texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)) {}

  // in
  vec2 gl_FragCoord;
  vec2 texCoord;

  void main()
  {
    float& W = _texelFetch(coverageTex, ivec2(gl_FragCoord));  //texture2D(coverageTex, vec2(gl_FragCoord)/view_wh);
    gl_FragColor = vec4(vec3(color), min(W, 1.0f)*color.a);
    // clear for any overlapping path
    W = 0.0f;
  }
};

} // end namespace CoverageFill


class CoverageSWRenderer : public SWRenderer
{
public:
  Buffer2D<float> coverageBuffer;
  CoverageAccum coverageAccum;

  CoverageSWRenderer() : coverageBuffer(new float[800*800], 800, 800), coverageAccum(coverageBuffer) { name = "SW Coverage"; }

  void drawGlyph(LBGlyph* g) override
  {
    ShaderPipeline<CoverageCalc::VS, CoverageCalc::FS, CoverageCalc::Varying, CoverageAccum> CoverageCalcPipeline(coverageAccum);
    CoverageCalc::verts = g->triangle_vertices;
    drawArrays(CoverageCalcPipeline, g->num_triangle_vertices);

    CoverageFill::verts = &quadVertices[0];
    CoverageFill::bbox_center = 0.5f*vec2(g->x_max + g->x_min, g->y_max + g->y_min);
    CoverageFill::bbox_wh = vec2(g->x_max - g->x_min, g->y_max - g->y_min);
    CoverageFill::coverageTex = coverageBuffer;
    ShaderPipeline<CoverageFill::VS, CoverageFill::FS, CoverageFill::Varying, BlendOver> CoverageFillPipeline(blendOver);
    drawArrays(CoverageFillPipeline, 6);
  }

  void beginFrame() override
  {
    coverageBuffer.fill(0);
    SWRenderer::beginFrame();
  }
};
