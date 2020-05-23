#include "swrender.h"


class MinDist {
public:
  Buffer2D<float>& distBuffer;
  MinDist(Buffer2D<float>& buff) : distBuffer(buff) {}

  // we expect FS to write distance to gl_FragColor.a
  void blend(int x, int y, const FSBase& fs)
  {
    float& dest = distBuffer.data[x + distBuffer.width*y];
    // note that we can accomplish this in OpenGL with glBlendEquation(GL_MIN)
    float src = fs.gl_FragColor.a;
    //dest = glm::min(dest, src);
    if(abs(src) < abs(dest))
       dest = src;
  }
};

// shader to generate a signed distance field - must be called for every pixel in SDF for every triangle
/*namespace SDFGenOld {

using namespace glm;

struct Varying : public VaryingBase {
  vec3 edge;
  vec3 edgescale;
};

const Vertex* verts;

const float ex = 12.0f; // extrusion amount, in pixels

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position((float)verts[idx].x, (float)verts[idx].y, (float)verts[idx].z), normal(verts[idx].triangle_normal),
    fringe((float)verts[idx].u, (float)verts[idx].v, (float)verts[idx].w), // gl_VertexID(idx),
    edge(vout->edge), edgescale(vout->edgescale) {}

  // in
  vec3 position;
  const vec2& normal;
  vec3 fringe;

  // out
  vec3& edge;
  vec3& edgescale;

  void main()
  {
    //float s2 = model[2][2];
    //edge = 0.5f*fringe*s2 + 1.0f*scale;

    edge = (0.5f/scale)*fringe;
    edgescale = abs(edge);
    // at most one of fringe.x,y,z should be < 0
    if(fringe.x < 0)
      edge.x -= ex*position.z;
    else if(fringe.x > 0)
      edge.x += ex;

    if(fringe.y < 0)
      edge.y -= ex*position.z;
    else if(fringe.y > 0)
      edge.y += ex;

    if(fringe.z < 0)
      edge.z -= ex*position.z;
    else if(fringe.z > 0)
      edge.z += ex;

    // expand triangle
    gl_Position = model * vec4(vec2(position) + ex*scale*normal, 0.0f, 1.0f);
  }
};

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    edge(baryinterp(vA.edge, vB.edge, vC.edge, t)),
    edgescale(baryinterp(vA.edgescale, vB.edgescale, vC.edgescale, t)) {}

  // in
  vec3 edge;
  vec3 edgescale;

  void main()
  {
    if(edgescale.x == 0 && edgescale.y == 0 && edgescale.z == 0) {
      gl_FragColor.a = -0.5f;
      return;
    }
    //vec3 vv = edgescale - edge;  //abs(edge)
    //float v = min(vv.x > 0 ? vv.x : FLT_MAX, min(vv.y > 0 ? vv.y : FLT_MAX, vv.z > 0 ? vv.z : FLT_MAX));
    vec3 vv = -edge + edgescale;
    float v = min(edgescale.x > 0 ? vv.x : FLT_MAX, min(edgescale.y > 0 ? vv.y : FLT_MAX, edgescale.z > 0 ? vv.z : FLT_MAX));
    gl_FragColor.a = -v;
  }
};

} // end namespace SDFGen
*/

namespace SDFGen {

using namespace glm;

//struct Vertex { vec2 position; vec2 normal; float fringe; };
struct Varying : public VaryingBase {
  vec2 va, vb, vc;
  vec3 edges;
};

const Vertex* verts;
bool drawTris;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx),
    position_in((float)verts[idx].x, (float)verts[idx].y), normal_in(verts[idx].edge_normal),
    edges_in((float)verts[idx].u, (float)verts[idx].v, (float)verts[idx].w), // gl_VertexID(idx),
    va_in(verts[idx].va), vb_in(verts[idx].vb), vc_in(verts[idx].vc),
    edges(vout->edges), va(vout->va), vb(vout->vb), vc(vout->vc) {}

  // in
  vec2 position_in;
  const vec2& normal_in;
  vec3 edges_in;
  const vec2& va_in;
  const vec2& vb_in;
  const vec2& vc_in;

  // out
  vec3& edges;

  vec2& va;
  vec2& vb;
  vec2& vc;

  void main()
  {
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
    //gl_FrontFacing(cross(vec2(vB.position) - vec2(vA.position), vec2(vC.position) - vec2(vA.position)) < 0),
    gl_FragCoord(t.x + 0.5f, t.y + 0.5f, 0.0f, 1.0f),
    edges(baryinterp(vA.edges, vB.edges, vC.edges, t)),
    va(baryinterp(vA.va, vB.va, vC.va, t)),
    vb(baryinterp(vA.vb, vB.vb, vC.vb, t)),
    vc(baryinterp(vA.vc, vB.vc, vC.vc, t)) {}

  //bool gl_FrontFacing;
  vec4 gl_FragCoord;

  // in
  vec3 edges;

  vec2 va;
  vec2 vb;
  vec2 vc;

  // distance from point `p` to line segment `a`-`b`
  float sdSegment(vec2 p, vec2 a, vec2 b)
  {
    vec2 ab = b - a;
    vec2 ap = p - a;
    float t = max(0.0f, min(1.0f, dot(ap, ab)/dot(ab, ab)));
    return distance(a + t*ab, p) * sign(ab.x*ap.y - ab.y*ap.x);
  }

  void main()
  {
    vec2 p(gl_FragCoord);
    float d_ab = sdSegment(p, va, vb);
    float d_bc = sdSegment(p, vb, vc);
    float d_ca = sdSegment(p, vc, va);

    float d = -1.0;  //-min(distance(p, va), min(distance(p, vb), distance(p, vc)));  // nearest point (always on path)
    if(edges.x)
      d = max(d, d_ab);
    if(edges.y)
      d = max(d, d_bc);
    if(edges.z)
      d = max(d, d_ca);
    if(abs(d) > 0.5f && (d_ab > 0 || d_bc > 0 || d_ca > 0))
      d = FLT_MAX; // exterior point far from edge

    if(drawTris)
      d = cos(3.14*d);

    gl_FragColor = vec4(0, 0, 0, d);


//    float a_ab = edges.x ? clamp(0.5f - sdSegment(p, va, vb), 0.0f, 1.0f) : 1.0f;
//    float a_bc = edges.y ? clamp(0.5f - sdSegment(p, vb, vc), 0.0f, 1.0f) : 1.0f;
//    float a_ca = edges.z ? clamp(0.5f - sdSegment(p, vc, va), 0.0f, 1.0f) : 1.0f;
//    float a_pt = 0.5f + min(distance(p, va), min(distance(p, vb), distance(p, vc)));
//    float a = min(min(a_ab, a_bc), min(a_ca, a_pt));

//    gl_FragColor = vec4(0, 0, 0, 0.5f - a);
  }
};

} // end namespace SDFGen


namespace SDFGen_LB {

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
    dFdx_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.x, t.db.x, t.dc.x)),
    dFdy_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.y, t.db.y, t.dc.y)) {}

  // in
  vec3 texCoord;
  vec3 dFdx_texCoord;
  vec3 dFdy_texCoord;

  void main()
  {
    vec2 p, px, py;
    float fx,fy,sd; //,alpha;

    p = vec2(texCoord);
    px = vec2(dFdx_texCoord);  // dFdx(p);
    py = vec2(dFdy_texCoord);  // dFdy(p);
    fx = 2.0f*p.x*px.x - px.y;
    fy = 2.0f*p.x*py.x - py.y;
    sd = (p.x*p.x  - p.y) / sqrt(fx*fx + fy*fy);
    gl_FragColor.a = sd * texCoord.z;
  }
};

} // end namespace SDFGen_LB


namespace SDF {

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
Buffer2D<float> distTex;

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FragCoord(t.x, t.y), texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)),
    dFdx_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.x, t.db.x, t.dc.x)),
    dFdy_texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t.da.y, t.db.y, t.dc.y)) {}

  // in
  ivec2 gl_FragCoord;
  vec2 texCoord;
  vec2 dFdx_texCoord;
  vec2 dFdy_texCoord;

  void main()
  {
    //anisotropicAA_deriv();
    //anisotropicAA(gl_FragCoord.x, gl_FragCoord.y);
    //distToAlpha(gl_FragCoord.x, gl_FragCoord.y);
    isotropicAA(gl_FragCoord.x, gl_FragCoord.y);
  }

  void anisotropicAA_deriv()
  {
    float D = texture2D(distTex, texCoord).r;
    //float aastep = length(vec2(dFdx(D), dFdy(D)));
    float dDdx = texture2D(distTex, texCoord + dFdx_texCoord).r - D;
    float dDdy =  texture2D(distTex, texCoord + dFdy_texCoord).r - D;
    // limit to 2 to handle FLT_MAX values outside of processed region
    // actually, change in D over a single pixel can't be more than sqrt(2)
    float aastep1 = min(2.0f, length(vec2(dDdx, dDdy)));
    float alpha = 1.0f - smoothstep(-aastep1, aastep1, D);
    gl_FragColor = vec4(vec3(color), color.a*alpha);
  }

  void anisotropicAA(int x, int y)
  {
    int dx = x+1 < distTex.width ? 1 : 0;
    int dy = y+1 < distTex.height ? 1 : 0;
    float D = distTex.data[x + y*distTex.width];
    float dDdx = distTex.data[x+dx + y*distTex.width] - D;
    float dDdy = distTex.data[x + (y+dy)*distTex.width] - D;
    float aastep1 = min(2.0f, length(vec2(dDdx, dDdy)));
    float alpha = 1.0f - smoothstep(-aastep1, aastep1, D);
    gl_FragColor = vec4(vec3(color), color.a*alpha);
  }

  //void distToAlpha(int x, int y)
  //{
  //  float D = distTex.data[x + y*distTex.width];
  //  gl_FragColor = vec4(vec3(color), D > 0 ? clamp((2.0f + SDFGenOld::ex - D)/20.0f, 0.0f, 1.0f) : 1.0f);
  //}

  void isotropicAA(int x, int y)
  {
    float& D = distTex.data[x + y*distTex.width];
    gl_FragColor = vec4(vec3(color), color.a*clamp(0.5f - D, 0.0f, 1.0f));

    // clear for next path - obviously can't do this here for anisotropic AA
    // to do this in OpenGL, we should be able to use image load/store (4.2 or extension)
    D = FLT_MAX;
  }
};

} // end namespace SDF


class SDFSWRenderer : public SWRenderer
{
public:
  Buffer2D<float> distBuffer;
  MinDist minDist;

  SDFSWRenderer() : distBuffer(new float[800*800], 800, 800), minDist(distBuffer) { name = "SW SDF 2-pass"; }

  // setup SDF renderer
  void drawGlyph(LBGlyph* g) override
  {
    // textures
    // Fill_FS::tex01 = &texture01[0];
    // vertices
    SDFGen::verts = g->triangle_vertices;

    ShaderPipeline<SDFGen::VS, SDFGen::FS, SDFGen::Varying, MinDist> SDFGenPipeline(minDist);
    drawArrays(SDFGenPipeline, g->num_triangle_vertices);

    SDF::verts = &quadVertices[0];
    SDF::bbox_center = 0.5f*vec2(g->x_max + g->x_min, g->y_max + g->y_min);
    SDF::bbox_wh = vec2(g->x_max - g->x_min, g->y_max - g->y_min);
    // distance field as texture
    SDF::distTex = distBuffer;
    ShaderPipeline<SDF::VS, SDF::FS, SDF::Varying, BlendOver> SDFPipeline(blendOver);
    drawArrays(SDFPipeline, 6);

    // inside curves
    /* SDFGen_LB::color = vec3(0, 0, 0);
    SDFGen_LB::verts = g->inside_curve_triangles;
    ShaderPipeline<SDFGen_LB::VS, SDFGen_LB::FS, SDFGen_LB::Varying, MinDist> SDFLBGen(minDist, ShaderInterface::SIGNED_DIST);
    drawArrays(SDFLBGen, g->num_inside_ctrl_points*3);

    // outside curves
    SDFGen_LB::verts = g->outside_curve_triangles;
    drawArrays(SDFLBGen, g->num_outside_ctrl_points*3); */
  }

  void beginFrame() override
  {
    distBuffer.fill(FLT_MAX);
    SWRenderer::beginFrame();
  }

//  void endFrame() override
//  {
//    SDF::verts = &quadVertices[0];
//    // distance field as texture
//    SDF::distTex = distBuffer;
//    ShaderPipeline<SDF::VS, SDF::FS, SDF::Varying, BlendOver> SDFPipeline(blendOver);
//    drawArrays(SDFPipeline, 6);

//    SWRenderer::endFrame();
//  }
};
