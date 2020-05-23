#include "swrender.h"

class WindingAccum {
public:
  Buffer2D<float>& windingBuffer;
  WindingAccum(Buffer2D<float>& buff) : windingBuffer(buff) {}

  // we expect FS to write distance to gl_FragColor.a
  void blend(int x, int y, const FSBase& fs)
  {
    float& dest = windingBuffer.data[x + windingBuffer.width*y];
    dest = dest + fs.gl_FragColor.a;
  }
};

// from sw_coverage.cpp
extern float coverageTri(vec2 vA, vec2 vB, vec2 vC, vec2 p);
extern float areaEdge(vec2 v0, vec2 v1);

/*
namespace WindingCalc {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 va, vb, vc;
};

const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
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
    gl_Position = model * vec4(position_in + 0.55f*scale*normal_in, 0.0f, 1.0f);
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

  // distance from point `p` to line segment `a`-`b`
  float sdSegment(const vec2& p, const vec2& a, const vec2& b)
  {
    vec2 ab = b - a;
    vec2 ap = p - a;
    float t = max(0.0f, min(1.0f, dot(ap, ab)/dot(ab, ab)));
    return distance(a + t*ab, p) * sign(ab.x*ap.y - ab.y*ap.x);
  }

  void main()
  {
    // for use with edge_normal (expand all edges)
    vec2 p = vec2(gl_FragCoord);
    float coverage = coverageTri(va, vb, vc, p);
    gl_FragColor = vec4(0, 0, 0, coverage);

    // for use with triangle_normal (expand only exterior edge)
//    vec2 p(gl_FragCoord);
//    float d = sdSegment(p, vb, vc);
//    //float w = gl_FrontFacing ? clamp(d + 0.5f, 0.0f, 1.0f) : clamp(d - 0.5f, -1.0f, 0.0f);
//    float w = clamp(d, -0.5f, 0.5f) + (gl_FrontFacing ? 0.5f : -0.5f);
//    gl_FragColor = vec4(0, 0, 0, w);

    // testing
//    gl_FragColor = gl_FrontFacing ? vec4(1, 0, 0, 0.5) : vec4(0, 0, 1, 0.5);
//    float w = clamp(0.5f - d, 0.0f, 1.0f);
//    if(abs(d) > 0.5 && gl_FrontFacing != (d > 0))
//      fprintf(stderr, "d: %f, front: %d\n", d, gl_FrontFacing ? 1 : -1);
//    gl_FragColor.a = gl_FrontFacing ? 1 : -1;
  }
};

} // end namespace WindingCalc
*/

// unlike areaEdge(), this assumes pixel center is (0,0), not (0.5, 0.5)
float areaEdge2(vec2 v0, vec2 v1)
{
  vec2 window = clamp(vec2(v0.x, v1.x), -0.5f, 0.5f);

  vec2 dv = v1 - v0;
  float slope = dv.y/dv.x;
  float midx = 0.5f*(window.x + window.y);
  float y = v0.y + (midx - v0.x)*slope;  // y value at middle of window

  float width = window.y - window.x;
  float dy = std::abs(slope*width);
  // credit for this to https://git.sr.ht/~eliasnaur/gio/tree/master/gpu/shaders/stencil.frag
  // if width == 1 (so midx == 0), the components of sides are: y crossing of right edge of frag, y crossing
  //  of left edge, x crossing of top edge, x crossing of bottom edge.  Since we only consider positive slope
  //  (note abs() above), there are five cases (below, bottom-right, left-right, left-top, above) - the area
  //  formula below reduces to these cases thanks to the clamping of the other values to 0 or 1.
  // I haven't thought carefully about the width < 1 case, but experimentally it matches areaEdge()
  vec4 sides = vec4(y + 0.5f*dy, y - 0.5f*dy, (0.5f - y)/dy, (-0.5f - y)/dy);  //ry, ly, tx, bx
  sides = clamp(sides + 0.5f, 0.0f, 1.0f);  // shift from -0.5..0.5 to 0..1 for area calc
  float area = 0.5f*(sides.z - sides.z*sides.y - 1.0f - sides.x + sides.x*sides.w);
  return width == 0.0f ? 0.0f : area * width;
}


namespace WindingTri {

using namespace glm;

struct Varying : public VaryingBase {};

const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position_in((float)verts[idx].x, (float)verts[idx].y), normal_in(verts[idx].edge_normal) {}

  // in
  vec2 position_in;
  const vec2& normal_in;

  void main()
  {
    gl_Position = model * vec4(position_in, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FrontFacing(cross(vec2(vB.position) - vec2(vA.position), vec2(vC.position) - vec2(vA.position)) > 0),
    gl_FragCoord(t.x + 0.5f, t.y + 0.5f, 0.0f, 1.0f) {}

  bool gl_FrontFacing;
  vec4 gl_FragCoord;

  void main()
  {
    gl_FragColor = gl_FrontFacing ? vec4(0, 0, 0, 1) : vec4(0, 0, 0, -1);
  }
};

} // end namespace WindingTri


// this takes trapezoid winding verts
namespace WindingEdge {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 va, vb, vc;
};

const TrapezoidVertex* verts;

vec2 origin;
float y_min;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    //position_in((float)verts[idx].x, (float)verts[idx].y),
    normal_in(verts[idx].triangle_normal),
    va_in(verts[idx].va), vb_in(verts[idx].vb),
    va(vout->va), vb(vout->vb) {}

  // in
  //vec2 position_in;
  const vec2& normal_in;
  const vec2& va_in;
  const vec2& vb_in;
  //const vec2& vc_in;

  // out
  vec2& va;
  vec2& vb;
  //vec2& vc;

  void main()
  {
    vec2 va_tf = vec2(model * vec4(va_in, 0.0f, 1.0f));
    vec2 vb_tf = vec2(model * vec4(vb_in, 0.0f, 1.0f));
    // screen coords for unexpanded vertices
    va = 0.5f*view_wh*(vec2(1.0f) + va_tf);
    vb = 0.5f*view_wh*(vec2(1.0f) + vb_tf);

    // normal.x > 0: we are right, otherwise left
    vec2 pos = (normal_in.x > 0) == (va_tf.x > vb_tf.x) ? va_tf : vb_tf;
    //vec2 pos = (normal_in.x > 0) ? vb_tf : va_tf;

    float y_max = max(va_tf.y, vb_tf.y);
    pos.y = normal_in.y > 0.0f ? y_max : y_min;

    //vec2 nba = normalize(vb_tf - va_tf);
    //vec2 expand = normal_in.x * nba + normal_in.y * vec2(-nba.y, nba.x);
    //PLATFORM_LOG("va, vb, normal, pos, expand: (%f, %f), (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n",
    //             va_tf.x, va_tf.y, vb_tf.x, vb_tf.y, normal_in.x, normal_in.y, pos.x, pos.y, expand.x, expand.y);

    // expand to get final position
    //gl_Position = vec4(pos + expand/view_wh, 0.0f, 1.0f);
    gl_Position = vec4(pos + 0.5f*(1.0f/400.0f)*normal_in, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FrontFacing(cross(vec2(vB.position) - vec2(vA.position), vec2(vC.position) - vec2(vA.position)) > 0),
    gl_FragCoord(t.x + 0.5f, t.y + 0.5f, 0.0f, 1.0f),
    va(baryinterp(vA.va, vB.va, vC.va, t)),
    vb(baryinterp(vA.vb, vB.vb, vC.vb, t)) {}

  bool gl_FrontFacing;
  vec4 gl_FragCoord;

  // in
  vec2 va;
  vec2 vb;

  void main()
  {
    vec2 p(gl_FragCoord);
    //float coverage = coverageTri(va, vb, origin, p);
    float coverage2 = areaEdge2(vb - p, va - p);

    p -= 0.5f;
    float coverage = va.x < vb.x ? areaEdge(va - p, vb - p) : -areaEdge(vb - p, va - p);

    if(std::abs(coverage2 - coverage) > 1E-5f)
      fprintf(stderr, "areaEdge: %f; gio: %f; diff: %f\n", coverage, coverage2, std::abs(coverage - coverage2)); //, to.x - from.x, to.y - from.y);

    gl_FragColor = vec4(0, 0, 0, coverage);
  }
};

} // end namespace WindingEdge

namespace WindingFill {

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
Buffer2D<float> windingTex;

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FragCoord(t.x + 0.5, t.y + 0.5), texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)) {}

  // in
  vec2 gl_FragCoord;
  vec2 texCoord;

  void main()
  {
    float& W = _texelFetch(windingTex, ivec2(gl_FragCoord));  //texture2D(windingTex, vec2(gl_FragCoord)/view_wh);
    //gl_FragColor = vec4(vec3(color), W == 0 ? 0 : color.a);
    gl_FragColor = vec4(vec3(color), min(abs(W), 1.0f)*color.a);
    // clear for any overlapping path
    W = 0.0f;
  }
};

} // end namespace WindingFill


class WindingSWRenderer : public SWRenderer
{
public:
  Buffer2D<float> windingBuffer;
  WindingAccum windingAccum;

  WindingSWRenderer() : windingBuffer(new float[800*800], 800, 800), windingAccum(windingBuffer) { name = "SW Winding"; }

  void drawGlyph(LBGlyph* g) override
  {
    //ShaderPipeline<WindingCalc::VS, WindingCalc::FS, WindingCalc::Varying, WindingAccum> WindingCalcPipeline(windingAccum);
    //WindingCalc::verts = g->winding_vertices;
    //drawArrays(WindingCalcPipeline, g->num_winding_vertices);

//    ShaderPipeline<WindingTri::VS, WindingTri::FS, WindingTri::Varying, WindingAccum> WindingTriPipeline(windingAccum);
//    for(int ii = 0; ii < g->num_contours; ++ii) {
//      WindingTri::verts = g->vertices + g->contours[ii];
//      int next_contour_idx = (ii + 1 < g->num_contours) ? g->contours[ii+1] : g->num_verts;
//      drawArrays(WindingTriPipeline, next_contour_idx - g->contours[ii], GL_TRIANGLE_FAN);
//    }

    ShaderPipeline<WindingEdge::VS, WindingEdge::FS, WindingEdge::Varying, WindingAccum> WindingEdgePipeline(windingAccum);

    // temporary hack - only works for single contour paths!
    vec4 origin_tf = UniformBase::model * vec4(g->vertices[0].x, g->vertices[0].y, 0, 1);
    WindingEdge::origin = 0.5f*UniformBase::view_wh*(vec2(1.0f) + vec2(origin_tf));

    vec4 c0 = UniformBase::model * vec4(g->x_min, g->y_min, 0, 1);
    vec4 c1 = UniformBase::model * vec4(g->x_max, g->y_min, 0, 1);
    vec4 c2 = UniformBase::model * vec4(g->x_min, g->y_max, 0, 1);
    vec4 c3 = UniformBase::model * vec4(g->x_max, g->y_max, 0, 1);
    WindingEdge::y_min = glm::min(glm::min(c0.y, c1.y), glm::min(c2.y, c3.y));

    WindingEdge::verts = g->trapezoid_vertices;
    drawArrays(WindingEdgePipeline, 6*g->num_verts, GL_TRIANGLES);

    WindingFill::verts = &quadVertices[0];
    WindingFill::bbox_center = 0.5f*vec2(g->x_max + g->x_min, g->y_max + g->y_min);
    WindingFill::bbox_wh = vec2(g->x_max - g->x_min, g->y_max - g->y_min);
    WindingFill::windingTex = windingBuffer;
    ShaderPipeline<WindingFill::VS, WindingFill::FS, WindingFill::Varying, BlendOver> WindingFillPipeline(blendOver);
    drawArrays(WindingFillPipeline, 6);
  }

  void beginFrame() override
  {
    windingBuffer.fill(0);
    SWRenderer::beginFrame();
  }
};
