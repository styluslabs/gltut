#include <vector>
#include "swrender.h"


/// Shader emulation

mat4 UniformBase::model;
float UniformBase::scale;
vec2 UniformBase::view_wh;
float UniformBase::aastep;
vec4 UniformBase::color;


static vec2 view_wh;
static vec2 view_xy;

void setViewport(int left, int top, int right, int bottom)
{
  view_xy = vec2(left, top);
  view_wh = vec2(right - left, bottom - top);
}

static glm::ivec2 debugPixel(-1,-1);
void setDebugPixel(int x, int y) { debugPixel = glm::ivec2(x, view_wh.y - y); }

// fragment coverage: if fragment is sufficient close to edge for possiblity of partial coverage (say any
//  barycentric coord, a, is within dadx+dady of 0), calculate intersection area
// can we be more precise in determining partial coverage? e.g., calculate barycentric coords of each
//  corner - partial coverage should not be possible if sign of coords is same for all corners

// summing partial coverage: setting alpha = coverage does not give correct result
// - it seems the only way to get correct result if we allow coverage calc for interior edges is to actually
//  keep track of the coverage mask/polygon.  One option would be to process each frag only once, i.e.,
//  process each polygon contributing to fragment in reverse z-order, stopping once coverage area reaches 1

// Sutherland-Hodgman polygon clipping algorithm - clipping polygon must be convex
//  ref: http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Lua
// Keep in mind that we need to handle the case of arbitrarily small triangles - we can't assume triangle
//  covers several fragments

bool halfPlane(const vec2& p, const vec2& cp1, const vec2& cp2)
{
  return (cp2.x-cp1.x)*(p.y-cp1.y) > (cp2.y-cp1.y)*(p.x-cp1.x);
}

vec2 intersection(const vec2& cp1, const vec2& cp2, const vec2& s, const vec2& e)
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

std::vector<vec2> clipPolygon(const std::vector<vec2>& subject, const std::vector<vec2>& clip)
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
// could summing the x coordinates cause loss of precision for polygons far from x axis?  perhaps calculation
//  using triangle fan (see nanovg) would be more robust?
float polygonArea(const std::vector<vec2>& points)
{
  float area = 0;
  for(size_t ii = 0, jj = points.size() - 1; ii < points.size(); jj = ii++)
    area += (points[jj].x + points[ii].x)*(points[jj].y - points[ii].y);
  return area/2;
}

// distance from point `pt` to line segment `start`-`end`
static float distToSegment(vec2 pt, vec2 start, vec2 end)
{
  const vec2 seg = end - start;
  if(start == end)
    return glm::distance(pt, start);
  // Consider the line extending the segment, parameterized as start + t*(end - start).
  // We find projection of pt onto this line and clamp t to [0,1] to limit to segment
  const float t = std::max(0.0f, std::min(1.0f, glm::dot(pt - start, end - start)/glm::dot(seg, seg)));
  const vec2 proj = start + t * (end - start);  // Projection falls on the segment
  return glm::distance(proj, pt);
}

// find closest point in triangle to given point p - from Real-time Collision Detection by Ericson
static float signedDistTriangle(const vec2& p, const vec2& a, const vec2& b, const vec2& c)
{
  // Check if P in vertex region outside A
  vec2 ab = b - a;
  vec2 ac = c - a;
  vec2 ap = p - a;
  float d1 = glm::dot(ab, ap);
  float d2 = glm::dot(ac, ap);
  if (d1 <= 0.0f && d2 <= 0.0f)
    return glm::distance(p, a); // barycentric coordinates (1,0,0)

  // Check if P in vertex region outside B
  vec2 bp = p - b;
  float d3 = glm::dot(ab, bp);
  float d4 = glm::dot(ac, bp);
  if (d3 >= 0.0f && d4 <= d3)
    return glm::distance(p, b); // barycentric coordinates (0,1,0)

  // Check if P in edge region of AB, if so return projection of P onto AB
  float vc = d1*d4 - d3*d2;
  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
    float v = d1 / (d1 - d3);
    return glm::distance(p, a + v * ab); // barycentric coordinates (1-v,v,0)
  }

  // Check if P in vertex region outside C
  vec2 cp = p - c;
  float d5 = glm::dot(ab, cp);
  float d6 = glm::dot(ac, cp);
  if (d6 >= 0.0f && d5 <= d6)
    return glm::distance(p, c); // barycentric coordinates (0,0,1)

  // Check if P in edge region of AC, if so return projection of P onto AC
  float vb = d5*d2 - d1*d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    float w = d2 / (d2 - d6);
    return glm::distance(p, a + w * ac); // barycentric coordinates (1-w,0,w)
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  float va = d3*d6 - d5*d4;
  vec2 bc = c - b;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return glm::distance(p, b + w * bc); // barycentric coordinates (0,1-w,w)
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  float denom = 1.0f/glm::abs(ab.x*ac.y - ab.y*ac.x);
  //float area_actual = glm::abs(ab.x*ac.y - ab.y*ac.x)/2.0f;
  //float area_bary = va*denom + vb*denom + vc*denom;
  return -glm::min(va*denom/glm::length(bc), glm::min(vb*denom/glm::length(ac), vc*denom/glm::length(ab)));

  /*float denom = 1.0f / (va + vb + vc);
  float v = vb * denom;
  float w = vc * denom;
  return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w */
}


vec2 toNDC(const vec4& v)
{
  return vec2(v)/v.w;
}

vec2 toScreen(const vec2& v)
{
  return (v + 1.0f)*view_wh/2.0f + view_xy;
}

float absmin(float a, float b)
{
  return glm::abs(a) < glm::abs(b) ? a : b;
}

float absmax(float a, float b)
{
  return glm::abs(a) > glm::abs(b) ? a : b;
}

void rasterizeTriangle(const ShaderInterface& prog, int vA, int vB, int vC)
{
  using glm::abs; using glm::min; using glm::max;

  // perspective division to give normalized device coordinates
  vec2 ndcA = toScreen(toNDC(prog.position(vA)));
  vec2 ndcB = toScreen(toNDC(prog.position(vB)));
  vec2 ndcC = toScreen(toNDC(prog.position(vC)));
  // triangle bounding box, clamped to viewport
  vec2 bbox_min = max(min(ndcA, min(ndcB, ndcC)), view_xy);
  vec2 bbox_max = min(max(ndcA, max(ndcB, ndcC)), view_xy + view_wh);
  if(bbox_min.x > bbox_max.x || bbox_min.y > bbox_max.y)
    return;
  // convert bbox from (-1..1, -1..1) to screen coords
  glm::ivec2 bbox_imin = glm::floor(bbox_min);
  glm::ivec2 bbox_imax = glm::ceil(bbox_max);

  // struct for barycentric coords and gradients
  FSInput t;
  // cartesian to barycentric transformation
  //mat2 Tinv = glm::inverse(mat2(ndcA - ndcC, ndcB - ndcC));
  //vec2 ab = Tinv*(((vec2(x,y) - view_xy)*2.0f/view_wh - 1.0f) - ndcC);  a = ab.x;  b = ab.y;
  // optimized version
  vec2 ndcCA = ndcA - ndcC;
  vec2 ndcCB = ndcB - ndcC;
  vec2 ndcCR = (vec2(bbox_imin) + 0.5f) - ndcC;
  float det = cross(ndcCB, ndcCA); // this is just 2*triangle area
  // calculate barycentric coords for bbox_min point
  float a_row = cross(ndcCB, ndcCR)/det;
  float b_row = cross(ndcCR, ndcCA)/det;
  // calculate step values
  t.da = vec2(-ndcCB.y/det, ndcCB.x/det);
  t.db = vec2(ndcCA.y/det, -ndcCA.x/det);
  t.dc = -(t.da + t.db);

  // iterate over bounding box pixels
  int ii, jj;
  for(t.y = bbox_imin.y, ii = 0; t.y < bbox_imax.y; ++t.y, ++ii) {
    t.a = a_row; t.b = b_row;
    for(t.x = bbox_imin.x, jj = 0; t.x < bbox_imax.x; ++t.x, ++jj) {
      // precision issues: ideally, we'd use ints for a,b,c - maybe introduce separate scale factors for each?
      //float a = a_row + dady*ii + dadx*jj, b = b_row + dbdy*ii + dbdx*jj;
      t.c = 1.0f - t.a - t.b;
      if(t.a > 0 && t.b > 0 && t.c > 0) {
        if(t.x == debugPixel.x && t.y == debugPixel.y)
          PLATFORM_LOG("Rendering pixel (%d,%d)\n", t.x, t.y);  // set breakpoint here

        prog.run_fs(vA, vB, vC, t);
      }
      t.a += t.da.x; t.b += t.db.x;  //c += t.dc.x;
    }
    a_row += t.da.y; b_row += t.db.y;
  }
}

/* Templated rasterizer

// some fixed point references:
// * http://stackoverflow.com/questions/79677/whats-the-best-way-to-do-fixed-point-math
// * http://wiki.yak.net/675/fixed.h
// * https://github.com/manuelbua/fixedpoint-math/blob/master/fixedpoint.h

template<unsigned int E>
struct Fixed_t
{
  typedef Fixed_t self;
  static constexpr int ONE = 1 << E;
  Fixed_t() : m(0) {}
  Fixed_t(int x) : m(x << E) {}
  Fixed_t(int x, int frac) : m((x << E) + frac) {}
  explicit Fixed_t(float d) : m(static_cast<int>(d * ONE + 0.5f)) {}
  explicit Fixed_t(double d) : m(static_cast<int>(d * ONE + 0.5)) {}

  explicit operator float() const { return float(m)/ONE; }
  explicit operator double() const { return double(m)/ONE; }
  explicit operator int() const { return m >> E; }

  self& operator+=(const self& x) { m += x.m; return *this; }
  self& operator-=(const self& x) { m -= x.m; return *this; }
  self& operator*=(const self& x) { m = (int)(((long long)m * x.m) >> E); return *this; }
  self& operator/=(const self& x) { m = (int)(((long long)m << E) / x.m); return *this; }
  self& neg() { m = -m; return *this; }
  self operator-() const { return self(*this).neg(); }
  bool operator!() const { return !m; }

  // using friend instead of member fns for binary ops allows left arg to be promoted via non-explicit constructor
  // Note: returning const value (to prevent e.g. (a + b) = ...) breaks C++11 move semantics
  friend self operator+(self x, const self& y) { return x += y; }
  friend self operator-(self x, const self& y) { return x -= y; }
  friend self operator*(self x, const self& y) { return x *= y; }
  friend self operator/(self x, const self& y) { return x /= y; }
  // comparison operators
  friend bool operator==(const self& x, const self& y) { return x.m == y.m; }
  friend bool operator!=(const self& x, const self& y) { return x.m != y.m; }
  friend bool operator>(const self& x, const self& y) { return x.m > y.m; }
  friend bool operator<(const self& x, const self& y) { return x.m < y.m; }
  friend bool operator>=(const self& x, const self& y) { return x.m >= y.m; }
  friend bool operator<=(const self& x, const self& y) { return x.m <= y.m; }

  int m;
};

// required for glm
namespace std {
  template<unsigned int E>
  struct numeric_limits< Fixed_t<E> >
  {
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = true;
    static constexpr bool is_iec559 = false;
  };
}

typedef Fixed_t<12> rscalar;
//typedef int rscalar;

// fixed point notes:
// - confirmed that incremental calc of a,b,c is always equal to multiplicative calc
//  - also, c = det-a-b is always equal to incremental calc (c += dc)
// - no missing pixels w/ doubles

// counterclockwise oriented A-B-C triangle

template<typename rscalar>
void rasterizeTriangleFixedPt(const ShaderInterface& prog, int vA, int vB, int vC)
{
  typedef glm::tvec2<rscalar, glm::highp> rvec2;
  using glm::abs; using glm::min; using glm::max; using glm::ceil; using glm::floor;

  // perspective division to give normalized device coordinates
  vec2 ndcA = toScreen(toNDC(prog.position(vA)));
  vec2 ndcB = toScreen(toNDC(prog.position(vB)));
  vec2 ndcC = toScreen(toNDC(prog.position(vC)));
  // triangle bounding box, clamped to viewport
  vec2 bbox_min = max(min(ndcA, min(ndcB, ndcC)), view_xy);
  vec2 bbox_max = min(max(ndcA, max(ndcB, ndcC)), view_xy + view_wh);
  if(bbox_min.x > bbox_max.x || bbox_min.y > bbox_max.y)
    continue;
  // convert bbox from (-1..1, -1..1) to screen coords
  glm::ivec2 bbox_imin = glm::floor(bbox_min);
  glm::ivec2 bbox_imax = glm::ceil(bbox_max);
  int w = bbox_imax.x - bbox_imin.x;

  // cartesian to barycentric transformation - all operations here must use rscalar!
  rvec2 ndcCA = rvec2(ndcA) - rvec2(ndcC);
  rvec2 ndcBC = rvec2(ndcC) - rvec2(ndcB);
  rvec2 ndcAB = rvec2(ndcB) - rvec2(ndcA);
  rvec2 ndcCR = (rvec2(bbox_imin) + rvec2(0.5f, 0.5f)) - rvec2(ndcC);
  rscalar det = cross(ndcCA, ndcBC); // this is just 2*triangle area
  float detf = float(det);

  // calculate barycentric coords for bbox_min point
  rscalar a_row = cross(ndcCR, ndcBC);
  rscalar b_row = cross(ndcCR, ndcCA);
  rscalar c_row = det - a_row - b_row; // cross(ndcAR, ndcAB);
  // calculate step values
  rvec2 da(ndcBC.y, -ndcBC.x);
  rvec2 db(ndcCA.y, -ndcCA.x);
  //rvec2 dc = -(da + db);  // (ndcAB.y, -ndcAB.x);

  // fill convention, edge included in triangle if on left or is horizontal top
  rscalar thresh_a = ndcBC.y < 0 || (ndcBC.y == 0 && ndcBC.x < 0) ? -rscalar(0, 1) : 0;
  rscalar thresh_b = ndcCA.y < 0 || (ndcCA.y == 0 && ndcCA.x < 0) ? -rscalar(0, 1) : 0;
  rscalar thresh_c = ndcAB.y < 0 || (ndcAB.y == 0 && ndcAB.x < 0) ? -rscalar(0, 1) : 0;

  // for determining range of x values within triangle for each scanline
  rscalar slope_a = da.y/da.x, slope_b = db.y/db.z, slope_c = dc.y/dc.x;
  rscalar ax = -a_row/da.x, bx = -b_row/db.x, cx = -c_row/dc.x;

  // struct for barycentric coords and gradients
  FSInput t;
  t.da = vec2(da)/detf;
  t.db = vec2(db)/detf;
  t.dc = -(t.da + t.db);

  // iterate over bounding box pixels
  for(int y = bbox_imin.y, ii = 0; y < bbox_imax.y; ++y, ++ii) {
    rscalar a = a_row, b = b_row;
    // calculate range of x values within triangle for this scanline
    // only perform this optimization if bbox width is above some minimum
    rscalar c = det - a - b;
    ax -= slope_a; bx -= slope_b; cx -= slope_c;
    // This is much, much, much, much, much uglier than I originally envisioned...
    int x0 = 0, x1 = w;
    // negative with pos d/dx
    if(a < 0 && ax > x0) x0 = ax;
    if(b < 0 && bx > x0) x0 = bx;
    if(c < 0 && cx > x0) x0 = cx;
    // positive with neg d/dx
    if(a > 0 && ax > 0 && ax < x1) x1 = ax+1;
    if(b > 0 && bx > 0 && bx < x1) x1 = bx+1;
    if(c > 0 && cx > 0 && cx < x1) x1 = cx+1;
    a += x0*da.x; b += x0*db.x;

    //for(int x = bbox_imin.x, jj = 0; x < bbox_imax.x; ++x, ++jj) {
    for(int x = bbox_imin.x + x0, jj = 0; x < bbox_imin.x + x1; ++x, ++jj) {
      c = det - a - b;
      if(a > thresh_a && b > thresh_b && c > thresh_c) {
        prog.run_fs(vA, vB, vC, t.update(x, y, float(a)/detf, float(b)/detf, float(c)/detf));
      }
      a += da.x; b += db.x;  //c += dc.x;
    }
    a_row += da.y; b_row += db.y;
  }
} */


void drawArrays(ShaderInterface& prog, int nverts, int geommode)
{
  prog.init(nverts);
  for(int ii = 0; ii < nverts; ++ii)
    prog.run_vs(ii);

  int incr = (geommode == GL_TRIANGLES) ? 3 : 1;
  for(int ii = 2; ii < nverts; ii += incr)
    rasterizeTriangle(prog, (geommode == GL_TRIANGLE_FAN ? 0 : ii-2), ii-1, ii);
}

// unlike glDrawElements, we take number of vertices as an arg and process them all up front
void drawElements(ShaderInterface& prog, int* indices, int nindices, int nverts, int geommode)
{
  prog.init(nverts);
  for(int ii = 0; ii < nverts; ++ii)
    prog.run_vs(ii);

  int incr = (geommode == GL_TRIANGLES) ? 3 : 1;
  for(int ii = 2; ii < nindices; ii += incr)
    rasterizeTriangle(prog, indices[(geommode == GL_TRIANGLE_FAN ? 0 : ii-2)], indices[ii-1], indices[ii]);
}

static const char* screenVert =
R"(
// below is a valid shader for #version 120
attribute vec2 position_in;
attribute vec2 texcoord_in;

varying vec2 texcoord;

void main()
{
  texcoord = texcoord_in;
  gl_Position = vec4(position_in, 0.0, 1.0);
}
)";

static const char* screenFrag =
R"(
// below is a valid shader for #version 120
varying vec2 texcoord;
uniform sampler2D texFramebuffer;

void main()
{
  gl_FragColor = texture2D(texFramebuffer, texcoord);
}
)";

unsigned int SWRenderer::clearColor = 0xFFFFFFFF;

SWRenderer::SWRenderer() : frameBuffer(new Color[800*800], 800, 800), blendOver(frameBuffer)
{
  spScreen.reset(new ShaderProgram(screenVert, screenFrag));
  glUniform1i(spScreen->uniform("texFramebuffer"), 0);
  glEnableVertexAttribArray(spScreen->attrib("position_in"));
  glEnableVertexAttribArray(spScreen->attrib("texcoord_in"));

  glGenTextures(1, &texScreen);
  glBindTexture(GL_TEXTURE_2D, texScreen);
  //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 800, 800, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //glGenVertexArrays(1, &vao);
}

SWRenderer::~SWRenderer()
{
  glDeleteTextures(1, &texScreen);
  delete[] frameBuffer.data;
  //glDeleteVertexArrays(1, &vao);
}

void SWRenderer::beginFrame()
{
  frameBuffer.fill<unsigned int>(clearColor);
}

void SWRenderer::endFrame()
{
  glUseProgram(spScreen->prog);
  //glBindVertexArray(vao);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texScreen);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 800, 800, 0, GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.data);
  glVertexAttribPointer(spScreen->attrib("position_in"), 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][0]);
  glVertexAttribPointer(spScreen->attrib("texcoord_in"), 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][2]);
  glDrawArrays(GL_TRIANGLES, 0, 6);
}
