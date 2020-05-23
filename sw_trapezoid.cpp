#include "swrender.h"

// identical to WindingAccum
class TrapezoidAccum {
public:
  Buffer2D<float>& windingBuffer;
  TrapezoidAccum(Buffer2D<float>& buff) : windingBuffer(buff) {}

  // we expect FS to write distance to gl_FragColor.a
  void blend(int x, int y, const FSBase& fs)
  {
    float& dest = windingBuffer.data[x + windingBuffer.width*y];
    dest = dest + fs.gl_FragColor.a;
  }
};

// from sw_coverage.cpp
extern float areaEdge(vec2 v0, vec2 v1);

namespace TrapezoidCalc {

using namespace glm;

struct Varying : public VaryingBase {
  vec2 va, vb, vc;
};

const TrapezoidVertex* verts;

// additional uniforms
float y_min;
bool uSide;

// shaders based on https://github.com/pcwalton/pathfinder/blob/master/shaders/gles2/stencil-aaa.*
class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    //position_in((float)verts[idx].x, (float)verts[idx].y),
    normal_in(verts[idx].triangle_normal),
    va_in(verts[idx].va), vb_in(verts[idx].vb), vc_in(verts[idx].vc),
    va(vout->va), vb(vout->vb), vc(vout->vc) {}

  // in
  //vec2 position_in;
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
    vec2 from = vec2(model * vec4(va_in, 0.0f, 1.0f));
    vec2 ctrl = vec2(model * vec4(vc_in, 0.0f, 1.0f));
    vec2 to = vec2(model * vec4(vb_in, 0.0f, 1.0f));

    vec2 v01 = ctrl - from, v12 = to - ctrl;
    float t = clamp(v01.x / (v01.x - v12.x), 0.0f, 1.0f);
    vec2 ctrl0 = mix(from, ctrl, t), ctrl1 = mix(ctrl, to, t);
    vec2 mid = mix(ctrl0, ctrl1, t);
    if(uSide) {
      from = mid;
      ctrl = ctrl1;
    }
    else {
      ctrl = ctrl0;
      to = mid;
    }
    // screen coords for unexpanded vertices
    va = 0.5f*view_wh*(vec2(1.0f) + from);
    vb = 0.5f*view_wh*(vec2(1.0f) + to);
    vc = 0.5f*view_wh*(vec2(1.0f) + ctrl);

    // if normal.x > 0 we are right, otherwise left
    vec2 pos = (normal_in.x > 0) == (from.x > to.x) ? from : to;

    float y_max = max(ctrl.y, max(from.y, to.y));
    pos.y = normal_in.y > 0 ? y_max : y_min;

    // expand to get final position
    gl_Position = vec4(pos + 0.5f*(1/400.0f)*normal_in, 0.0f, 1.0f);

    // zero area?
    if(abs(from.x - to.x) < 0.00001f)
      gl_Position = vec4(-2.0f, -2.0f, 0.0, 1.0);    // move off-screen

    // expand triangle by half a pixel so all partially covered pixels are included
    //gl_Position = model * vec4(position_in + 0.5f*scale*normal_in, 0.0f, 1.0f);
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
    p -= 0.5f;
    vec2 from = va - p, ctrl = vc - p, to = vb - p;

    // Determine winding, and sort into a consistent order so we only need to find one root below.
    bool winding = from.x < to.x;
    vec2 left = winding ? from : to, right = winding ? to : from;
    vec2 v0 = ctrl - left, v1 = right - ctrl;

    // find curve parameter `t` for x = <center of pixel>
    //vec2 window = clamp(vec2(from.x, to.x), -0.5f, 0.5f);
    vec2 window = clamp(vec2(from.x, to.x), 0.0f, 1.0f);
    float offset = mix(window.x, window.y, 0.5f) - left.x;
    float t = offset / (v0.x + sqrt(v1.x * offset - v0.x * (offset - v0.x)));

    // Compute position and derivative to form a line approximation.
    float y = mix(mix(left.y, ctrl.y, t), mix(ctrl.y, right.y, t), t);
    float d = mix(v0.y, v1.y, t) / mix(v0.x, v1.x, t);

    // Look up area under that line, and scale horizontally to the window size.
    float dX = window.x - window.y;

    vec2 r0(0, y - d*dX/2);
    vec2 r1(abs(dX), y + d*dX/2);
    float coverage = sign(dX) * areaEdge(r0, r1);

    float width = window.y - window.x;
    float dy = abs(d*width);
    y -= 0.5f;
    vec4 sides = vec4(dy*+0.5f + y, dy*-0.5f + y, (+0.5f-y)/dy, (-0.5f-y)/dy);
    sides = clamp(sides+0.5f, 0.0f, 1.0f);
    float area = 0.5f*(sides.z - sides.z*sides.y - 1.0f - sides.x+sides.x*sides.w);
    area *= width;  //winding ? -width : width;
    if(std::abs(area - coverage) > 1E-6f)
     fprintf(stderr, "areaEdge: %f; gio: %f; dx,dy: %f, %f\n", coverage, area, to.x - from.x, to.y - from.y);

    //float coverage = (va.x < vb.x ? areaEdge(va, vb) : -areaEdge(vb, va));
    gl_FragColor = vec4(0, 0, 0, coverage);


    //vec2 p(gl_FragCoord);
    //p -= 0.5f;
    //va -= p;  vb -= p;
    //float coverage = (va.x < vb.x ? areaEdge(va, vb) : -areaEdge(vb, va));
    //gl_FragColor = vec4(0, 0, 0, coverage);
  }
};

} // end namespace TrapezoidCalc

namespace TrapezoidFill {

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

} // end namespace TrapezoidFill

// this renderer computes fractional winding number for trapezoids between edge and bottom of glyph bbox
// - 1 areaEdge() call per pixel per edge vs. 3 for WindingSWRenderer, but two triangles per edge vs. 1 and
//  possibly more overdraw (thus less than 3x reduction in total areaEdge calls!)

class TrapezoidSWRenderer : public SWRenderer
{
public:
  Buffer2D<float> windingBuffer;
  TrapezoidAccum trapezoidAccum;

  TrapezoidSWRenderer() : windingBuffer(new float[800*800], 800, 800), trapezoidAccum(windingBuffer) { name = "SW Trapezoid"; }

  void drawGlyph(LBGlyph* g) override
  {
    ShaderPipeline<TrapezoidCalc::VS, TrapezoidCalc::FS, TrapezoidCalc::Varying, TrapezoidAccum> TrapezoidCalcPipeline(trapezoidAccum);
    TrapezoidCalc::verts = g->trapezoid_vertices;

    // bottom of trapezoids in transformed coords
    vec4 c0 = UniformBase::model * vec4(g->x_min, g->y_min, 0, 1);
    vec4 c1 = UniformBase::model * vec4(g->x_max, g->y_min, 0, 1);
    vec4 c2 = UniformBase::model * vec4(g->x_min, g->y_max, 0, 1);
    vec4 c3 = UniformBase::model * vec4(g->x_max, g->y_max, 0, 1);
    TrapezoidCalc::y_min = glm::min(glm::min(c0.y, c1.y), glm::min(c2.y, c3.y));

    TrapezoidCalc::uSide = false;
    drawArrays(TrapezoidCalcPipeline, g->num_trapezoid_vertices, GL_TRIANGLES);
    TrapezoidCalc::uSide = true;
    drawArrays(TrapezoidCalcPipeline, g->num_trapezoid_vertices, GL_TRIANGLES);

    TrapezoidFill::verts = &quadVertices[0];
    TrapezoidFill::bbox_center = 0.5f*vec2(g->x_max + g->x_min, g->y_max + g->y_min);
    TrapezoidFill::bbox_wh = vec2(g->x_max - g->x_min, g->y_max - g->y_min);
    TrapezoidFill::windingTex = windingBuffer;
    ShaderPipeline<TrapezoidFill::VS, TrapezoidFill::FS, TrapezoidFill::Varying, BlendOver> TrapezoidFillPipeline(blendOver);
    drawArrays(TrapezoidFillPipeline, 6);
  }

  void beginFrame() override
  {
    windingBuffer.fill(0);
    SWRenderer::beginFrame();
  }
};
