#include <stdio.h>

#include "lbfont.h"

// using GLU tessellator
#include <GL/glu.h>
#ifdef _MSC_VER
//#include <windows.h>
//#define GLU_TESS_CB void(CALLBACK*)()
#define GLU_TESS_CALL __stdcall
#define GLU_TESS_CB void(GLU_TESS_CALL*)()
#else
#define GLU_TESS_CALL
#define GLU_TESS_CB _GLUfuncptr
#endif

#define STB_TRUETYPE_IMPLEMENTATION  //  force following include to generate implementation
#include "stb_truetype.h"


int freadall(const char* filename, char** bufferout)
{
  FILE* f = fopen(filename, "rb");
  if(!f)
    return 0;
  // obtain file size
  fseek(f, 0, SEEK_END);
  int size = ftell(f);
  rewind(f);
  char* buf = (char*)malloc(size+1);
  if(!buf)
    return 0;
  *bufferout = buf;
  // copy the file into the buffer:
  int bytesread = fread(buf, 1, size, f);
  fclose(f);
  buf[bytesread] = '\0';
  return bytesread;
}

// nothing to do really...always going to get only GL_TRIANGLES because edgeFlag callback is set
void GLU_TESS_CALL tessBeginCB(GLenum which, LBGlyph* g) { }
void GLU_TESS_CALL tessEndCB(LBGlyph* g) { }
void GLU_TESS_CALL tessEdgeFlagCB(GLboolean flag, LBGlyph* g) { }

static int tessErrors = 0;
// hopefully will never get called, only error that is given by GLU
// could be that GLU_TESS_COMBINE callback is needed, i dont think that
// should actually be the case with the font outlines though
void GLU_TESS_CALL tessErrorCB(GLenum errorCode, LBGlyph*)
{
  const GLubyte *errorStr;
  errorStr = gluErrorString(errorCode);
  PLATFORM_LOG("Tesselator Error: %s\n", errorStr);
  ++tessErrors;
}

// add the index of the vertex to the indices array for drawing the triangles later
void GLU_TESS_CALL tessVertexCB(Vertex* vertex, LBGlyph* g)
{
  GLuint index = vertex - g->vertices;
  g->triangle_indices[g->num_triangle_indices++] = index;
}

// calculate vertex extrusion vector, given two edge normals
// assumes FS will multiply by 0.5px at current scale
static vec2 extrusionDir(const vec2& n1, const vec2& n2)
{
  // old approach - always a normalized vector
  //return glm::normalize(n1 + n2);
  // new approach - return proper length to give 0.5px expansion of edges
  // note that as vertex angle -> 0, the length of this vector -> infinity
  return (n1 + n2)/(1 + glm::dot(n1, n2));
}

// remember: tesslator doesn't know anything about lineto, curveto, etc. - it just considers the
//  passed vertices as forming a polygon
static void finishContour(LBGlyph* g, int v0, int vi)
{
  vec2 n1 = glm::normalize(vec2(-(g->vertices[v0].y - g->vertices[vi].y), g->vertices[v0].x - g->vertices[vi].x));
  g->vertices[vi].path_normal = extrusionDir(n1, g->vertices[vi].path_normal);
  // fix up next/prev vertex info
  g->vertices[v0].u = vi;  // wrap prev vertex index for first vertex
  g->vertices[vi].w = v0;  // wrap next vertex index for last vertex
}

static GLdouble prev_vertex_x, prev_vertex_y;
static void addVertex(LBGlyph* g, int vi, GLdouble x, GLdouble y, bool first)
{
  g->vertices[vi].x = x;
  g->vertices[vi].y = y;
  g->vertices[vi].z = 1.0;

  vec2 n1 = glm::normalize(vec2(-(y - prev_vertex_y), x - prev_vertex_x));
  g->vertices[vi].path_normal = n1;
  if(!first)
    g->vertices[vi-1].path_normal = extrusionDir(g->vertices[vi-1].path_normal, n1);
  prev_vertex_x = x; prev_vertex_y = y;

  // temporarily store index information in Vertex for use by buildTriangleData
  g->vertices[vi].u = vi - 1;
  g->vertices[vi].v = vi;
  g->vertices[vi].w = vi + 1;
}

// tesselate a single glyph from the font (whichever one is loaded into teh glyph slot)
int tesselateOutline(LBGlyph* g, stbtt_vertex* points, int n_points, float EM_per_unit)
{
  // Create and configure a new tessellation object
  GLUtesselator* tess = gluNewTess();
  gluTessCallback(tess, GLU_TESS_ERROR_DATA,     (GLU_TESS_CB)tessErrorCB);
  gluTessCallback(tess, GLU_TESS_BEGIN_DATA,     (GLU_TESS_CB)tessBeginCB);
  gluTessCallback(tess, GLU_TESS_END_DATA,       (GLU_TESS_CB)tessEndCB);
  gluTessCallback(tess, GLU_TESS_EDGE_FLAG_DATA, (GLU_TESS_CB)tessEdgeFlagCB);
  gluTessCallback(tess, GLU_TESS_VERTEX_DATA,    (GLU_TESS_CB)tessVertexCB);
  tessErrors = 0;

  float sx = EM_per_unit;
  float sy = EM_per_unit;
  int vi = 0, v0 = 0, contour_index = 0;
  gluTessBeginPolygon(tess, g);
  gluTessBeginContour(tess);
  for(int i = 0; i < n_points; i++) {
    // we expect first point to always be a moveto command
    if(points[i].type == STBTT_vmove) {
      if(i > 0) {
        gluTessEndContour(tess);
        gluTessBeginContour(tess);
        finishContour(g, v0, vi-1);
      }
      v0 = vi;
      g->contours[contour_index] = vi;
      contour_index++;
      prev_vertex_x = points[i].x * sx;
      prev_vertex_y = points[i].y * sy;
      continue;
    }

    if(points[i].type == STBTT_vcurve) {
      // control points are only part of the tesselation if on the inside of the shape
      // check if its left or right (in/out) between last and next point
      int ax, ay, bx, by, x, y;
      addVertex(g, vi, points[i].cx * sx, points[i].cy * sy, vi == v0);
      x = points[i].cx; y = points[i].cy;
      ax = points[i-1].x; ay = points[i-1].y;
      bx = points[i].x; by = points[i].y;
      if((bx-ax)*(y-ay) - (by-ay)*(x-ax) < 0){
        gluTessVertex(tess, &g->vertices[vi].x, &g->vertices[vi].x);
        g->inside_ctrl_indices[g->num_inside_ctrl_points++] = vi;
      }
      else
        g->outside_ctrl_indices[g->num_outside_ctrl_points++] = vi;
      // use invalid index value to prevent antialiasing of triangles for interior control point
      // only necessary for interior points, but no harm in doing for outer points too!
      g->vertices[vi].v = -1;
      vi++;
    }
    addVertex(g, vi, points[i].x * sx, points[i].y * sy, vi == v0);
    // first arg is pointer to vertex data, 2nd is userdata passed to GLU_TESS_VERTEX_DATA callback
    gluTessVertex(tess, &g->vertices[vi].x, &g->vertices[vi].x);
    vi++;
  }
  finishContour(g, v0, vi-1);

  gluTessEndContour(tess);
  gluTessEndPolygon(tess);
  gluDeleteTess(tess);
  return tessErrors;
}

static bool areAdjacent(const Vertex* a, const Vertex* b)
{
  // u = prev, v = curr, w = next
  return a->v >= 0 && b->v >= 0 && (a->w == b->v || b->w == a->v);
}

// clockwise normal
static vec2 normal(const vec2& v)
{
  return glm::normalize(vec2(-v.y, v.x));
}

void buildWindingData(LBGlyph* g)
{
  //Vertex origin;
  //float eps = 0.01f;
  //origin.x = 0.5*(g->x_min + g->x_max);  //- eps;
  //origin.y = 0.5*(g->y_min + g->y_max);  //- eps;

  int vidx = 0;
  for(int ii = 0; ii < g->num_contours; ++ii) {
    int next_contour_idx = (ii + 1 < g->num_contours) ? g->contours[ii+1] : g->num_verts;
    int origin_idx = g->contours[ii];
    for(int jj = origin_idx + 1; jj < next_contour_idx; ++jj) {
      int next_vertex_idx = (jj + 1 < next_contour_idx) ? jj+1 : g->contours[ii];
      g->winding_vertices[vidx++] = g->vertices[origin_idx];
      g->winding_vertices[vidx++] = g->vertices[jj];
      g->winding_vertices[vidx++] = g->vertices[next_vertex_idx];
    }
  }

  for(int ii = 0; ii < vidx; ii += 3) {
    Vertex* a = &g->winding_vertices[ii];
    Vertex* b = &g->winding_vertices[ii+1];
    Vertex* c = &g->winding_vertices[ii+2];
    vec2 pta(a->x, a->y);
    vec2 ptb(b->x, b->y);
    vec2 ptc(c->x, c->y);

    vec2 ab = ptb - pta;
    vec2 ac = ptc - pta;
    vec2 cb = ptb - ptc;
    //bool bc_adj = areAdjacent(b, c);
    bool frontfacing = cross(ptb - pta, ptc - pta) > 0;
    vec2 ncb = normal(frontfacing ? cb : -cb);

    a->triangle_normal = vec2(0, 0);
    b->triangle_normal = ab/glm::dot(ab, ncb);
    c->triangle_normal = ac/glm::dot(ac, ncb);

    // we'll use edge normal to expand triangle in all directions (incl. interior edges)
    vec2 nab = normal(frontfacing ? -ab : ab);
    vec2 nac = normal(frontfacing ? ac : -ac);
    a->edge_normal = extrusionDir(nab, nac);
    b->edge_normal = extrusionDir(nab, ncb);
    c->edge_normal = extrusionDir(nac, ncb);

    a->va = pta; a->vb = ptb; a->vc = ptc;
    b->va = pta; b->vb = ptb; b->vc = ptc;
    c->va = pta; c->vb = ptb; c->vc = ptc;
  }
  g->num_winding_vertices = vidx;
}


void buildTrapezoidData(LBGlyph* g, stbtt_vertex* points, int n_points, float EM_per_unit)
{
  TrapezoidVertex v0, v1, v2, v3;
  // v0 - v1 - v2 - v3 always clockwise from top left
  v0.triangle_normal = vec2(-1, 1);
  v1.triangle_normal = vec2(1, 1);
  // we should actually use y_min of each subpath
  //v2.y = g->y_min;
  //v3.y = g->y_min;
  v2.triangle_normal = vec2(1, -1);
  v3.triangle_normal = vec2(-1, -1);

  float sx = EM_per_unit;
  float sy = EM_per_unit;
  int vidx = 0;
  vec2 pta, ptb;  //contour_start;
  for(int i = 0; i < n_points; i++) {
    // we expect first point to always be a moveto command, and expect all contours to be closed
    if(points[i].type == STBTT_vmove) {
      ptb = vec2(points[i].x * sx, points[i].y * sy);;
      continue;
    }

    pta = ptb;
    ptb = vec2(points[i].x * sx, points[i].y * sy);
    // ptc is control point
    vec2 ptc;
    if(points[i].type == STBTT_vcurve)
      ptc = vec2(points[i].cx * sx, points[i].cy * sy);
    else
      ptc = 0.5f*(pta + ptb);

    v0.va = pta; v1.va = pta; v2.va = pta; v3.va = pta;
    v0.vb = ptb; v1.vb = ptb; v2.vb = ptb; v3.vb = ptb;
    v0.vc = ptc; v1.vc = ptc; v2.vc = ptc; v3.vc = ptc;

    // ideally, we'd use glMultiDrawElements or primitive restart to draw each trapezoid as a triangle strip
    //  at least, we should use glDrawElements (w/ GL_TRIANGLES)
    g->trapezoid_vertices[vidx++] = v0;
    g->trapezoid_vertices[vidx++] = v1;
    g->trapezoid_vertices[vidx++] = v2;
    g->trapezoid_vertices[vidx++] = v0;
    g->trapezoid_vertices[vidx++] = v2;
    g->trapezoid_vertices[vidx++] = v3;
  }
  g->num_trapezoid_vertices = vidx;
}

void buildTriangleData(LBGlyph* g)
{
  // because fringe value for a given vertex can be different for the different triangles sharing the
  //  vertex, we must have separate vertices for each triangle and use glDrawArrays instead of glDrawElements
  for(int ii = 0; ii < g->num_triangle_indices; ++ii)
    g->triangle_vertices[ii] = g->vertices[g->triangle_indices[ii]];

  g->num_triangle_vertices = g->num_triangle_indices;

  // now, we need to find exterior triangle edges (edges where vertices are consecutive)
  // we'll ignore the single triangle case; in all other cases, either one or two edges will be exterior
  for(int ii = 0; ii < g->num_triangle_indices; ii += 3) {
    Vertex* a = &g->triangle_vertices[ii];
    Vertex* b = &g->triangle_vertices[ii+1];
    Vertex* c = &g->triangle_vertices[ii+2];
    vec2 pta(a->x, a->y);
    vec2 ptb(b->x, b->y);
    vec2 ptc(c->x, c->y);

    vec2 ab = ptb - pta;
    vec2 ac = ptc - pta;
    vec2 cb = ptb - ptc;
    // need to calculate and save these before .v and .w are overwritten
    bool ab_adj = areAdjacent(a, b);
    bool bc_adj = areAdjacent(b, c);
    bool ca_adj = areAdjacent(c, a);

    // calculate extrusion vectors to expand each edge of triangle for SDF generation
    vec2 nab = normal(ab);
    vec2 nac = normal(-ac);
    vec2 ncb = normal(-cb);

    if(ab_adj && ca_adj)
      a->triangle_normal = extrusionDir(nab, nac);
    else if(ab_adj)
      a->triangle_normal = ac/glm::dot(ac, nab);
    else if(ca_adj)
      a->triangle_normal = ab/glm::dot(ab, nac);
    else
      a->triangle_normal = vec2(0, 0);

    if(ab_adj && bc_adj)
      b->triangle_normal = extrusionDir(nab, ncb);
    else if(ab_adj)
      b->triangle_normal = cb/glm::dot(cb, nab);
    else if(bc_adj)
      b->triangle_normal = ab/glm::dot(ab, ncb);
    else
      b->triangle_normal = vec2(0, 0);

    if(ca_adj && bc_adj)
      c->triangle_normal = extrusionDir(nac, ncb);
    else if(ca_adj)
      c->triangle_normal = cb/glm::dot(cb, nac);
    else if(bc_adj)
      c->triangle_normal = ac/glm::dot(ac, ncb);
    else
      c->triangle_normal = vec2(0, 0);

    // we'll use edge normal to expand triangle in all directions (incl. interior edges)
    a->edge_normal = extrusionDir(nab, nac);
    b->edge_normal = extrusionDir(nab, ncb);
    c->edge_normal = extrusionDir(nac, ncb);

    // use z component of vertex to store factor for determining effect of extrusion on height
    a->z = -glm::dot(ncb, a->triangle_normal);
    b->z = -glm::dot(nac, b->triangle_normal);
    c->z = -glm::dot(nab, c->triangle_normal);

    a->va = pta; a->vb = ptb; a->vc = ptc;
    b->va = pta; b->vb = ptb; b->vc = ptc;
    c->va = pta; c->vb = ptb; c->vc = ptc;

    a->u = ab_adj; a->v = bc_adj; a->w = ca_adj;
    b->u = ab_adj; b->v = bc_adj; b->w = ca_adj;
    c->u = ab_adj; c->v = bc_adj; c->w = ca_adj;


//    float twicearea = glm::abs(ab.x*ac.y - ab.y*ac.x);
//    float hc = twicearea/glm::length(ab);
//    a->u = ab_adj ? hc : -hc; b->u = ab_adj ? hc : -hc; c->u = -hc;
//    float ha = twicearea/glm::length(cb);
//    a->v = -ha; b->v = bc_adj ? ha : -ha; c->v = bc_adj ? ha : -ha;
//    float hb = twicearea/glm::length(ac);
//    a->w = ca_adj ? hb : -hb; b->w = -hb; c->w = ca_adj ? hb : -hb;

//    // default to no AA
//    a->u = 0; b->u = 0; c->u = 0;
//    a->v = 0; b->v = 0; c->v = 0;
//    a->w = 0; b->w = 0; c->w = 0;
//    //continue;  // disable triangle AA
//    // calculate area of triangle to aid point-line distance calc (2D cross-product of two edge vectors)
//    float twicearea = glm::abs(ab.x*ac.y - ab.y*ac.x);
//    // we need to write this value to an attribute array; shader will use it for antialiasing
//    if(ab_adj) {
//      // segment ab is on exterior; calculate distance from line ab to point c
//      float h = twicearea/glm::length(ab);
//      a->u = h; b->u = h; c->u = -h;
//    }
//    if(bc_adj) {
//      float h = twicearea/glm::length(cb);
//      a->v = -h; b->v = h; c->v = h;
//    }
//    if(ca_adj) {
//      float h = twicearea/glm::length(ac);
//      a->w = h; b->w = -h; c->w = h;
//    }
  }
}

void buildCurveData(LBGlyph* g)
{
  int contour_index=0, inside_index=0, outside_index=0;
  for(int i = 0; i < g->num_verts; ++i) {
    int before = i-1;
    int after = i+1;
    // wrap 'before' to last vertex of contour
    if(i == g->contours[contour_index])
      before = contour_index+1 < g->num_contours ? g->contours[contour_index+1] - 1 : g->num_verts - 1;
    // wrap 'after' to first vertex of contour, and advance to next countour if necessary
    if(after == g->num_verts)
      after = g->contours[contour_index];
    else if(contour_index+1 < g->num_contours && after == g->contours[contour_index+1]) {
      after = g->contours[contour_index];
      ++contour_index;
    }

    if(inside_index < g->num_inside_ctrl_points && i == (int)g->inside_ctrl_indices[inside_index]) {
      Vertex a,b,c;
      a = g->vertices[after];  // point after
      b = g->vertices[before]; // point before
      c = g->vertices[i];      // control point

      // calculate extrusion vectors to expand each edge of triangle for SDF generation
      vec2 pta(a.x, a.y);
      vec2 ptb(b.x, b.y);
      vec2 ptc(c.x, c.y);
      vec2 nab = normal(ptb - pta);
      vec2 nac = normal(pta - ptc);
      vec2 ncb = normal(ptc - ptb);
      a.triangle_normal = extrusionDir(nab, nac);
      b.triangle_normal = extrusionDir(nab, ncb);
      c.triangle_normal = extrusionDir(nac, ncb);

      // compute texture coords as values to quadratic equation for solving bezier
      // misusing red channel to mark inside vs outside for now
      a.u = 0.0;  a.v = 0.0;  a.w = -1.0;
      b.u = 1.0;  b.v = 1.0;  b.w = -1.0;
      c.u = 0.5;  c.v = 0.0;  c.w = -1.0;

      g->inside_curve_triangles[3*inside_index+0] = a;
      g->inside_curve_triangles[3*inside_index+1] = b;
      g->inside_curve_triangles[3*inside_index+2] = c;
      inside_index++;
    }

    if(outside_index < g->num_outside_ctrl_points && i == (int)g->outside_ctrl_indices[outside_index]) {
      Vertex a,b,c;
      a = g->vertices[after];  // point after
      b = g->vertices[before]; // point before
      c = g->vertices[i];      // control point

      // calculate extrusion vectors to expand each edge of triangle for SDF generation
      vec2 pta(a.x, a.y);
      vec2 ptb(b.x, b.y);
      vec2 ptc(c.x, c.y);
      vec2 nab = normal(ptb - pta);
      vec2 nac = normal(pta - ptc);
      vec2 ncb = normal(ptc - ptb);
      a.triangle_normal = extrusionDir(nab, nac);
      b.triangle_normal = extrusionDir(nab, ncb);
      c.triangle_normal = extrusionDir(nac, ncb);

      // compute texture coords as values to quadratic equation for solving bezier
      a.u = 0.0;  a.v = 0.0;  a.w = 1.0;
      b.u = 1.0;  b.v = 1.0;  b.w = 1.0;
      c.u = 0.5;  c.v = 0.0;  c.w = 1.0;

      g->outside_curve_triangles[3*outside_index+0] = a;
      g->outside_curve_triangles[3*outside_index+1] = b;
      g->outside_curve_triangles[3*outside_index+2] = c;
      outside_index++;
    }
  }
}

// for miter join
/* static int writeStrokeVertex(const Vertex& in, Vertex* out)
{
  out[0] = in;
  out[0].u = 1;
  out[1] = in;
  out[1].u = -1;
  return 2;
} */

// I don't think we can use miter join for inside angle because it will give incorrect result for very short segments
// for bevel and round join
static int writeStrokeVertex(const Vertex& prev, const Vertex& curr, const Vertex& next, Vertex* out)
{
  out[0] = curr;
  out[0].u = 1;  out[0].v = 1;
  out[0].edge_normal = glm::normalize(vec2(-(prev.y - curr.y), prev.x - curr.x));
  out[1] = curr;
  out[1].u = -1;  out[1].v = 1;
  out[1].edge_normal = out[0].edge_normal;

  out[2] = curr;
  out[2].u = -1;  out[2].v = -1;
  out[2].edge_normal = glm::normalize(-vec2(-(next.y - curr.y), next.x - curr.x));
  out[3] = curr;
  out[3].u = 1;  out[3].v = -1;
  out[3].edge_normal = out[2].edge_normal;

  return 4;
}

void buildStrokeData(LBGlyph* g)
{
  g->stroke_vertices = (Vertex*)calloc(4*(g->num_verts + g->num_contours), sizeof(Vertex));
  Vertex* v0 = g->vertices;  // this will be pointer to first vertex of contour
  int vidx = 0;
  for(int ii = 0; ii < g->num_contours; ++ii) {
    int nelem = ((ii+1 < g->num_contours) ? g->contours[ii+1] : g->num_verts) - g->contours[ii];
    vidx += writeStrokeVertex(v0[nelem-1], v0[0], v0[1], &g->stroke_vertices[vidx]);
    for(int jj = 1; jj < nelem - 1; ++jj)
      vidx += writeStrokeVertex(v0[jj-1], v0[jj], v0[jj+1], &g->stroke_vertices[vidx]);
    vidx += writeStrokeVertex(v0[nelem-2], v0[nelem-1], v0[0], &g->stroke_vertices[vidx]);
    // repeat first vertex to close path
    vidx += writeStrokeVertex(v0[nelem-1], v0[0], v0[1], &g->stroke_vertices[vidx]);
    // advance to next countour
    v0 += nelem;
  }
}

void initLBGlyphPts(LBGlyph* g, float scale, int n_points, stbtt_vertex* points)
{
  int n_ctrl_pts = 0;
  float x0 = 0, y0 = 0;
  for(int i = 0; i < n_points; i++) {
    // build bounding rect
    float x = scale*points[i].x, y = scale*points[i].y;
    g->x_min = i == 0 ? x : glm::min(x, g->x_min);  g->x_max = i == 0 ? x : glm::max(x, g->x_max);
    g->y_min = i == 0 ? y : glm::min(y, g->y_min);  g->y_max = i == 0 ? y : glm::max(y, g->y_max);

    if(points[i].type == STBTT_vmove)
      ++g->num_contours;
    else if(i > 0)
      g->area += (x + x0)*(y0 - y)/2;
    x0 = x;  y0 = y;

    if(points[i].type == STBTT_vcurve) {
      ++n_ctrl_pts;
      float cx = scale*points[i].cx, cy = scale*points[i].cy;
      g->x_min = glm::min(cx, g->x_min);  g->x_max = glm::max(cx, g->x_max);
      g->y_min = glm::min(cy, g->y_min);  g->y_max = glm::max(cy, g->y_max);
    }
  }
  // subtract num_contours since paths are explicity closed (last point == first point), but we don't want
  //  to pass duplicate vertices to tesselator
  g->num_verts = n_points + n_ctrl_pts - g->num_contours;

  // contours starts the starting index in vertices for each contour
  g->contours = (int*)calloc(g->num_contours, sizeof(int));
  g->vertices = (Vertex*)calloc(g->num_verts, sizeof(Vertex));
  g->triangle_indices = (GLuint*)calloc(3*g->num_verts, sizeof(GLuint)); // can be more, up to 3*#verts
  g->inside_ctrl_indices = (GLuint*)calloc(n_ctrl_pts, sizeof(GLuint)); // worst case...all are inside
  g->outside_ctrl_indices = (GLuint*)calloc(n_ctrl_pts, sizeof(GLuint)); // worst case...all are outside

  // parse outline/tesselate the glyph
  if(tesselateOutline(g, points, n_points, scale) > 0)
    PLATFORM_LOG("Errors encountered while tesslating '%c'\n", g->c);

  // generate triangles
  // worst case is every triangle having two exterior edges (assuming none have 3)
  g->triangle_vertices = (Vertex*)calloc(2 * g->num_triangle_indices, sizeof(Vertex));
  buildTriangleData(g);

  // allocate memory for curve vertices and compute border triangles/texture coords
  g->inside_curve_triangles  = (Vertex*)calloc(g->num_inside_ctrl_points * 3, sizeof(Vertex));
  g->outside_curve_triangles = (Vertex*)calloc(100, sizeof(Vertex));
  buildCurveData(g);

  g->winding_vertices = (Vertex*)calloc(3*g->num_verts + 1, sizeof(Vertex));
  buildWindingData(g);

  // two triangles per edge (= 2 tris per vertex)
  g->trapezoid_vertices = (TrapezoidVertex*)calloc(6*n_points, sizeof(TrapezoidVertex));
  buildTrapezoidData(g, points, n_points, scale);

  // must duplicate vertices for stroking
  buildStrokeData(g);
}

LBGlyph* initLBGlyph(LBFont* font, char c)
{
  int cidx = stbtt_FindGlyphIndex(font->fontinfo, (int)c);
  int advanceWidth;
  int leftSideBearing;
  stbtt_GetGlyphHMetrics(font->fontinfo, cidx, &advanceWidth, &leftSideBearing);
  // get outline
  stbtt_vertex* points;
  int n_points = stbtt_GetGlyphShape(font->fontinfo, cidx, &points);

  LBGlyph* g = (LBGlyph*)calloc(1, sizeof(LBGlyph));
  g->c = c;
  g->advance = (float)advanceWidth*(float)font->EM_per_unit;
  if(n_points > 0)
    initLBGlyphPts(g, font->EM_per_unit, n_points, points);
  stbtt_FreeShape(font->fontinfo, points);
  return g;
}

void freeLBGlyph(LBGlyph* g)
{
  if(!g) return;
  if(g->num_verts > 0) {
    free(g->contours);
    free(g->vertices);
    free(g->triangle_indices);
    free(g->inside_ctrl_indices);
    free(g->outside_ctrl_indices);
    free(g->inside_curve_triangles);
    free(g->outside_curve_triangles);
    free(g->triangle_vertices);
    free(g->stroke_vertices);
    free(g->winding_vertices);
    free(g->trapezoid_vertices);
  }
  free(g);
}

void drawGlyphData(LBGlyph* g)
{
  glEnableClientState(GL_VERTEX_ARRAY);
  glPointSize(5.0);

  // draw a blue dot for all vertices in the glyph
  glColor3f(0,0,1);
  glVertexPointer(3, LB_REAL, sizeof(Vertex), &g->vertices[0].x);
  glDrawArrays(GL_POINTS, 0, g->num_verts);

  // draw a smaller red dot for all inside bezier control points
  glColor3f(1,0,0);
  glDrawElements(GL_POINTS, g->num_inside_ctrl_points, GL_UNSIGNED_INT, g->inside_ctrl_indices);

  // draw a smaller green dot for all outside bezier control points
  glColor3f(0,1,0);
  glDrawElements(GL_POINTS, g->num_outside_ctrl_points, GL_UNSIGNED_INT, g->outside_ctrl_indices);

  glDisableClientState(GL_VERTEX_ARRAY);
}

LBFont* createLBFont(const char* filename)
{
  LBFont* self = (LBFont*)malloc(sizeof(LBFont));
  self->ttfbuffer = NULL;
  if(freadall(filename, &self->ttfbuffer) == 0) {
    free(self);
    return NULL;
  }
  self->fontinfo = (stbtt_fontinfo*)malloc(sizeof(stbtt_fontinfo));
  stbtt_InitFont(self->fontinfo, (unsigned char*)self->ttfbuffer, 0);
  self->EM_per_unit = stbtt_ScaleForMappingEmToPixels(self->fontinfo, 1);

  // load glyphs
  self->num_glyphs = 128;
  for(int i = 0; i < self->num_glyphs; i++)
    self->glyphs[i] = i < 32 ? NULL : initLBGlyph(self, i);

  return self;
}

void freeLBFont(LBFont* font)
{
  for(int ii = 0; ii < font->num_glyphs; ++ii)
    freeLBGlyph(font->glyphs[ii]);
  free(font->fontinfo);
  free(font->ttfbuffer);
  free(font);
}
