#pragma once

#include <algorithm>  // required for std::min,max on Windows
#include <glm/glm.hpp>

#include "glpp.h"

// typedef as an alternative to importing glm namespace - also allows us to use doubles instead of floats later
typedef glm::vec2 vec2;
typedef glm::vec3 vec3;
typedef glm::vec4 vec4;
typedef glm::mat2 mat2;
typedef glm::mat4 mat4;
typedef glm::ivec2 ivec2;

// glm only provides cross for vec3
template<typename T>
T cross(const glm::tvec2<T, glm::highp>& a, const glm::tvec2<T, glm::highp>& b)
{
  return a.x*b.y - a.y*b.x;
}

// GLU tessellator expects doubles
#define LB_REAL GL_DOUBLE
typedef struct {
  GLdouble x, y, z;
  GLdouble u, v, w;
  vec2 va, vb, vc;
  vec2 path_normal;
  vec2 triangle_normal;
  vec2 edge_normal;
  //GLdouble s, t;
} Vertex;

// minimal Vertex for performance testing
typedef struct {
  vec2 va, vb, vc;
  vec2 triangle_normal;
} TrapezoidVertex;

typedef struct {
  char c;
  float advance;
  float x_min, y_min, x_max, y_max;
  float area;
  int num_verts;
  int num_contours;
  int num_triangle_indices;
  int num_triangle_vertices;
  int num_inside_ctrl_points;
  int num_outside_ctrl_points;
  int num_winding_vertices;
  int num_trapezoid_vertices;
  int* contours;
  Vertex* vertices;
  Vertex* triangle_vertices;
  Vertex* inside_curve_triangles;
  Vertex* outside_curve_triangles;
  Vertex* stroke_vertices;
  Vertex* winding_vertices;
  TrapezoidVertex* trapezoid_vertices;
  GLuint* triangle_indices;
  GLuint* inside_ctrl_indices;
  GLuint* outside_ctrl_indices;
} LBGlyph;

class Renderer
{
public:
  //Renderer(const char* name) : strName(name) {}
  virtual ~Renderer() {}
  virtual void drawGlyph(LBGlyph* g) = 0;
  virtual void beginFrame() {}
  virtual void endFrame() {}
  virtual bool isSWRenderer() const { return false; }

  const char* name = "Abstract Renderer";
};

// we are free to add additional uniforms to individual classes
class UniformBase {
public:
  // VS uniforms
  static mat4 model;
  static float scale; // units per pixel
  static vec2 view_wh;
  // FS uniforms
  static float aastep;
  static vec4 color;
};
