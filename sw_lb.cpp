#include <vector>
#include <memory>

#include "swrender.h"
#include "glpp.h"


// note: namespace must be named - anonymous namespace will allow enclosed "using namespace ..." to break out
namespace Fill {

using namespace glm;

//struct Vertex { vec2 position; vec2 normal; float fringe; };
struct Varying : public VaryingBase {
  vec2 va, vb, vc;
  vec3 edges;
  //vec3 edge;
  //vec3 edgescale;
};

const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position_in((float)verts[idx].x, (float)verts[idx].y), normal_in(verts[idx].triangle_normal),
    edges_in((float)verts[idx].u, (float)verts[idx].v, (float)verts[idx].w), // gl_VertexID(idx),
    va_in(verts[idx].va), vb_in(verts[idx].vb), vc_in(verts[idx].vc),
    edges(vout->edges), va(vout->va), vb(vout->vb), vc(vout->vc) {}
    //edge(vout->edge), edgescale(vout->edgescale), origvertex(vout->origvertex) {}

  // in
  vec2 position_in;
  const vec2& normal_in;
  //vec3 fringe;
  vec3 edges_in;
  const vec2& va_in;
  const vec2& vb_in;
  const vec2& vc_in;

  // out
  //vec3& edge;
  //vec3& edgescale;
  vec3& edges;

  vec2& va;
  vec2& vb;
  vec2& vc;

  void main()
  {
    //float s2 = 0.5f/scale;
    //edge = fringe*s2;
    //edgescale = abs(edge);
    //if(edge.x < 0) edge.x += -0.5;
    //if(edge.y < 0) edge.y += -0.5;
    //if(edge.z < 0) edge.z += -0.5;


    // screen coords for unexpanded vertices
    va = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(va_in, 0.0f, 1.0f)));
    vb = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vb_in, 0.0f, 1.0f)));
    vc = 0.5f*view_wh*(vec2(1.0f)+vec2(model * vec4(vc_in, 0.0f, 1.0f)));

    edges = edges_in;

    // expand triangle by half a pixel so all partially covered pixels are included
    gl_Position = model * vec4(position_in + 0.5f*scale*normal_in, 0.0f, 1.0f);
  }
};


// I think gl_FrontFacing is correct for GL_CCW - note that y-axis is down, so z-axis is into screen for right-handed coords
class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FrontFacing(cross(vec2(vB.position) - vec2(vA.position), vec2(vC.position) - vec2(vA.position)) < 0),
    gl_FragCoord(t.x + 0.5f, t.y + 0.5f, 0.0f, 1.0f),
    edges(baryinterp(vA.edges, vB.edges, vC.edges, t)),
    va(baryinterp(vA.va, vB.va, vC.va, t)),
    vb(baryinterp(vA.vb, vB.vb, vC.vb, t)),
    vc(baryinterp(vA.vc, vB.vc, vC.vc, t)) {}
//    edge(baryinterp(vA.edge, vB.edge, vC.edge, t)),
//    edgescale(baryinterp(vA.edgescale, vB.edgescale, vC.edgescale, t)) {}

  bool gl_FrontFacing;
  vec4 gl_FragCoord;

  // in
  //vec3 edge;
  //vec3 edgescale;
  vec3 edges;

  vec2 va;
  vec2 vb;
  vec2 vc;

  // http://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
  // negative inside
  float sdTriangle(vec2 p, vec2 p0, vec2 p1, vec2 p2)
  {
    vec2 e0 = p1-p0, e1 = p2-p1, e2 = p0-p2;
    vec2 v0 = p -p0, v1 = p -p1, v2 = p -p2;

    vec2 pq0 = v0 - e0*clamp(dot(v0,e0)/dot(e0,e0), 0.0f, 1.0f);
    vec2 pq1 = v1 - e1*clamp(dot(v1,e1)/dot(e1,e1), 0.0f, 1.0f);
    vec2 pq2 = v2 - e2*clamp(dot(v2,e2)/dot(e2,e2), 0.0f, 1.0f);

    float s = sign( e0.x*e2.y - e0.y*e2.x );
    vec2 d = min(min(vec2(dot(pq0,pq0), s*(v0.x*e0.y-v0.y*e0.x)),
                     vec2(dot(pq1,pq1), s*(v1.x*e1.y-v1.y*e1.x))),
                     vec2(dot(pq2,pq2), s*(v2.x*e2.y-v2.y*e2.x)));

    return -sqrt(d.x)*sign(d.y);
  }

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
    vec2 p(gl_FragCoord);
    //float d = sdTriangle(vec2(gl_FragCoord), va, vb, vc);
    float a_ab = edges.x ? clamp(0.5f - sdSegment(p, va, vb), 0.0f, 1.0f) : 1.0f;
    float a_bc = edges.y ? clamp(0.5f - sdSegment(p, vb, vc), 0.0f, 1.0f) : 1.0f;
    float a_ca = edges.z ? clamp(0.5f - sdSegment(p, vc, va), 0.0f, 1.0f) : 1.0f;
    float a_pt = 0.5f + min(distance(p, va), min(distance(p, vb), distance(p, vc)));
    float a = min(min(a_ab, a_bc), min(a_ca, a_pt));
    //float a = clamp(0.5f - d, 0.0f, 1.0f);
    /*
    // remember that value of edge here is interpolated from 3 vertices defining our triangle
    // note that we could shift this calculation to vs, letting interpolation generate alpha directly
    //float a = edgescale > 0 ? gl_FragCoverage : (gl_FragCoverage >= 0.5 ? 1.0f : 0.0f);
    //gl_FragColor = vec4(color, a);

    vec3 vv = edgescale - edge;  //abs(edge)

    float v = min(vv.x > 0 ? vv.x : 1.0f, min(vv.y > 0 ? vv.y : 1.0f, vv.z > 0 ? vv.z : 1.0f));
//    float a = v > 0 ? smoothstep(0.0f, aastep, v) : 1.0f;
    float a = clamp(v, 0.0f, 1.0f);

    if(!gl_FrontFacing)
      a = 0;

    if(v < 1.0f)
      a = clamp((v < 0.5f ? 0.5f - dist : 0.5f + dist), 0.0f, 1.0f);

//    float v = min(vv.x, min(vv.y, vv.z));
//    float a = smoothstep(0.0f, aastep, v);
    */

    gl_FragColor = vec4(vec3(color), a*color.a);
    // this is handy for debugging
//    if((edgescale.x > 0 && edgescale.x < 0.5*aastep) || (edgescale.y > 0 && edgescale.y < 0.5*aastep) || (edgescale.z > 0 && edgescale.z < 0.5*aastep))
//      gl_FragColor = vec4(0, 0, 0, 1);

//    else if(abs(edge.z) < edgescale.z)
//      gl_FragColor = mix(vec4(1.0, 0.0, 0.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), 0.5 + 0.5*edge.z/edgescale.z);
//    else
//      gl_FragColor = mix(vec4(1.0, 1.0, 0.0, 1.0), vec4(0.0, 0.5, 0.0, 1.0), 0.5 + 0.5*edge.y/edgescale.y);
  }
};

} // end namespace Fill

namespace LB {

using namespace glm;

struct Varying : public VaryingBase {
  vec3 texCoord;
};

// putting uniforms outside shader classes mirrors the fact that GL VS and FS uniform names share the same namespace
const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position((float)verts[idx].x, (float)verts[idx].y, (float)verts[idx].z),
    texcoord((float)verts[idx].u, (float)verts[idx].v, (float)verts[idx].w),
    texCoord(vout->texCoord) {}

  // in
  vec3 position;
  vec3 texcoord;

  // out
  vec3& texCoord;

  void main()
  {
    texCoord = texcoord;
    gl_Position = model * vec4(position.x, position.y, 0.0, 1.0);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)),
    // Other options for implementing dFdx:
    // - run FS in 2x2 blocks, setup values in constructor, allow run() to access neighbors' values
    dFdx_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.x, t.db.x, t.dc.x)),
    dFdy_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.y, t.db.y, t.dc.y)) {}

  // in
  vec3 texCoord;
  vec3 dFdx_texCoord;
  vec3 dFdy_texCoord;

  void main()
  {
    vec2 p, px, py;
    float fx,fy,sd,alpha;

    p = vec2(texCoord);
    px = vec2(dFdx_texCoord);  // dFdx(p);
    py = vec2(dFdy_texCoord);  // dFdy(p);
    fx = 2.0f*p.x*px.x - px.y;
    fy = 2.0f*p.x*py.x - py.y;
    sd = (p.x*p.x  - p.y) / sqrt(fx*fx + fy*fy);
    alpha = 0.5f - sd * texCoord.z;
    alpha = clamp(alpha, 0.0f, 1.0f);
    gl_FragColor = vec4(vec3(color), alpha*color.a);
  }
};

} // end namespace LB

class LBSWRenderer : public SWRenderer
{
public:
  LBSWRenderer() { name = "SW SDF 1-pass"; }

  void drawGlyph(LBGlyph* g) override
  {
    // textures
    // Fill_FS::tex01 = &texture01[0];
    // vertices
    Fill::verts = g->triangle_vertices;

    ShaderPipeline<Fill::VS, Fill::FS, Fill::Varying, BlendOver> FillPipeline(blendOver);
    drawArrays(FillPipeline, g->num_triangle_vertices);

    // inside curves
    LB::verts = g->inside_curve_triangles;
    ShaderPipeline<LB::VS, LB::FS, LB::Varying, BlendOver> LBPipeline(blendOver);
    drawArrays(LBPipeline, g->num_inside_ctrl_points*3);

    // outside curves
    LB::verts = g->outside_curve_triangles;
    drawArrays(LBPipeline, g->num_outside_ctrl_points*3);
  }
};
