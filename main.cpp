#include <iostream>
#include <vector>
#include <random>
#include <functional>

#define GLEW_STATIC
#include <GL/glew.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

// Eigen is a popular alternative to GLM for more advanced linear algebra
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "glpp.h"

#include "tesselator.h"


const char* extrudeVS =
R"(#version 120

attribute vec2 position;
attribute vec2 normal;
attribute float miter;
//uniform mat4 proj;
uniform mat4 model;
//uniform mat4 view;
//uniform float thickness;
varying float edge;
varying float edgescale;

void main()
{
  edge = miter;
  edgescale = abs(miter);
  vec2 pointPos = position.xy + vec2(normal * miter);  // * thickness/2.0
  gl_Position = /* proj * view * */ model * vec4(pointPos, 0.0, 1.0);
  //gl_PointSize = 1.0;
}
)";

const char* extrudeFS =
R"(#version 120

uniform vec3 color;
varying float edge;
varying float edgescale;

void main()
{
  // remember that value of edge here is interpolated from 3 vertices defining our triangle
  // note that we could shift this calculation to vs, letting interpolation generate alpha directly
  float v = edgescale - abs(edge);
  v = smoothstep(0.0, 0.004, v); // v*inner);
  gl_FragColor = vec4(color, v);  //mix(vec4(color, 1.0), vec4(0.0), v);
}
)";


const char* extrude2VS =
R"(#version 120

attribute vec2 position;
attribute vec2 normal;
attribute float miter;

uniform mat4 model;

varying float edge;

void main()
{
  edge = sign(miter);
  vec2 pointPos = position.xy + vec2(normal * miter);
  gl_Position = /* proj * view * */ model * vec4(pointPos, 0.0, 1.0);
}
)";

const char* extrude2FS =
R"(#version 120

uniform vec3 color;
varying float edge;

void main()
{
  float dist = abs(edge);
  float afwidth = 1.4 * length(vec2(dFdx(dist), dFdy(dist)));
  float aa = smoothstep(0.0, afwidth, 1.0 - dist);
  gl_FragColor = vec4(color, aa);
}
)";


const char* fillVS =
R"(#version 120

attribute vec2 position;
//attribute vec2 normal;
attribute float fringe;

uniform mat4 model;

varying float edge;
varying float edgescale;

void main()
{
  float scale = model[2][2];
  edge = fringe*scale;
  edgescale = abs(fringe*scale);
  gl_Position = /* proj * view * */ model * vec4(position.xy, 0.0, 1.0);
}
)";

const char* fillFS =
R"(#version 120

uniform vec3 color;
varying float edge;
varying float edgescale;

void main()
{
  // remember that value of edge here is interpolated from 3 vertices defining our triangle
  // note that we could shift this calculation to vs, letting interpolation generate alpha directly
  float v = edgescale - abs(edge);
  v = edgescale > 0 ? smoothstep(0.0, 0.004, v) : 1.0;
  gl_FragColor = vec4(color, v);
  // this is handy for debugging
  //gl_FragColor = mix(vec4(0.0, 1.0, 0.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0), 0.5 + 0.5*edge/edgescale);
}
)";


// fill
int main(int argc, char* argv[])
{
  SDL_Init(SDL_INIT_VIDEO);

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);

  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow("OpenGL", 100, 100, 800, 600, SDL_WINDOW_OPENGL);
  SDL_GLContext context = SDL_GL_CreateContext(window);

  glewExperimental = GL_TRUE;
  glewInit();

  // libtess2 tesselator
  TESStesselator* tess = 0;
  tess = tessNewTess(NULL); //&ma);
  if(!tess)
    return -1;

  //VertexBuffer* vbScreen = new VertexBuffer(quadVertices, sizeof(quadVertices), QUAD_VERT_SIZE);
  //VertexArray* vaScreen = new VertexArray;
  std::vector<glm::vec2> points;
  std::vector<glm::vec2> normals;
  std::vector<glm::vec2> vertex;
  std::vector<GLfloat> fringe;
  std::vector<int> triangles;  // triplets of indices into points defining triangles
  auto addPoint = [&](int xi, int yi) {
    float x = (xi - 400.0f)/400.0f;
    float y = -(yi - 300.0f)/300.0f;
    int N = points.size();
    float prevx = N > 0 ? points[N-1].x : 0.0f;
    float prevy = N > 0 ? points[N-1].y : 0.0f;
    // require minimum spacing between points
    if(std::abs(x - prevx) + std::abs(y - prevy) < 0.02)
      return;
    points.push_back(glm::vec2(x, y));
    if(N > 0) {
      glm::vec2 n1 = glm::normalize(glm::vec2(-(y - prevy), x - prevx));
      if(N > 1) {
        normals[N-1] = glm::normalize(normals[N-1] + n1);
        // it seems contour must be explicitly closed for libtess2
        points.push_back(points.front());
        tessAddContour(tess, 2, &points[0], sizeof(glm::vec2), points.size());
        points.pop_back();
        if(!tessTesselate(tess, TESS_WINDING_POSITIVE, TESS_POLYGONS, 3, 2, 0))
          printf("tessTesselate ERROR!");
        //const float* verts = tessGetVertices(tess);
        const int* vinds = tessGetVertexIndices(tess);
        const int* elems = tessGetElements(tess);
        //const int nverts = tessGetVertexCount(tess);
        const int nelems = tessGetElementCount(tess);
        // convert indices into libtess2's vertex array into indices into our vertex array (points)
        int nverts = 3*nelems;
        triangles.resize(nverts);
        for(int ii = 0; ii < nverts; ++ii)
          triangles[ii] = vinds[elems[ii]];

        // because fringe value for a given vertex can be different for the different triangles sharing the
        //  vertex, we must have separate vertex input values for each triangle and use glDrawArrays instead of glDrawElements
        vertex.resize(nverts);
        for(int ii = 0; ii < nverts; ++ii)
          vertex[ii] = points[triangles[ii]];

        // now, we need to find exterior triangle edges (edges where vertices are consecutive)
        // we'll ignore the single triangle case; in all other cases, either one or two edges will be exterior

        // dealing with 2 edge case:
        // - split into two triangles?
        // - provide two AA params per vertex and mix or take min/max
        fringe.resize(nverts);
        for(int ii = 0; ii < nverts; ii += 3) {
          int a = triangles[ii];
          int b = triangles[ii+1];
          int c = triangles[ii+2];
          glm::vec2 ab = points[b] - points[a];
          glm::vec2 ac = points[c] - points[a];
          glm::vec2 cb = points[b] - points[c];
          // default to no AA
          fringe[ii] = 0; fringe[ii+1] = 0; fringe[ii+2] = 0;
          // calculate area of triangle to aid point-line distance calc (2D cross-product of two edge vectors)
          float twicearea = glm::abs(ab.x*ac.y - ab.y*ac.x);
          // we need to write this value to an attribute array; shader will use it for antialiasing
          if(std::abs(a - b) == 1 || std::abs(a - b) == nverts - 1) {
            // segment ab is on exterior; calculate distance from line ab to point c
            float h = twicearea/glm::length(ab);
            fringe[ii] = h; fringe[ii+1] = h; fringe[ii+2] = -h;
          }
          if(std::abs(b - c) == 1 || std::abs(b - c) == nverts - 1) {
            float h = twicearea/glm::length(cb);
            fringe[ii] = -h; fringe[ii+1] = h; fringe[ii+2] = h;
          }
          if(std::abs(c - a) == 1 || std::abs(c - a) == nverts - 1) {
            float h = twicearea/glm::length(ac);
            fringe[ii] = h; fringe[ii+1] = -h; fringe[ii+2] = h;
          }
          //else
          //  printf("Triangle has no exterior edges!\n");
        }
      }
      else
        normals.push_back(n1);
      normals.push_back(n1);
    }
  };

  ShaderProgram* spScreen = new ShaderProgram(fillVS, fillFS);
  //vaScreen->bindAttrib(spScreen, vbScreen, "position", 2, 0);
  //vaScreen->bindAttrib(spScreen, vbScreen, "texcoord", 2, 2);

  GLint attrPos = glGetAttribLocation(spScreen->prog, "position");
  glEnableVertexAttribArray(attrPos);
  // If an attribute or uniform is not used in shader, it may be optimized away and trying to access it will fail
  GLint attrFringe = glGetAttribLocation(spScreen->prog, "fringe");
  glEnableVertexAttribArray(attrFringe);

  // VS uniforms
  GLuint uniModel = spScreen->uniform("model");
  // FS uniforms
  glUniform3f(spScreen->uniform("color"), 0.0, 0.0, 1.0);

  //vaScreen->bind();
  spScreen->use();

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  checkGLError();

  glm::mat4 basexform;
  glm::mat4 xform;
  int initX, initY;
  SDL_Event event;
  bool drawTris = false;
  while(1) {
    if(SDL_PollEvent(&event)) {
      if(event.type == SDL_QUIT)
        break;
      else if(event.type == SDL_KEYUP) {
        if(event.key.keysym.sym == SDLK_ESCAPE)
          break;
        else if(event.key.keysym.sym == SDLK_F1)
          drawTris = !drawTris;
      }
      else if(event.type == SDL_MOUSEBUTTONDOWN) {
        if(event.button.button == SDL_BUTTON_LEFT) {
          points.clear();
          normals.clear();
          xform = glm::mat4();
          basexform = glm::mat4();
          addPoint(event.button.x, event.button.y);
        }
        else if(event.button.button == SDL_BUTTON_RIGHT) {
          initX = event.button.x;
          initY = event.button.y;
          basexform = xform;
        }
      }
      else if(event.type == SDL_MOUSEMOTION) {
        if(event.motion.state & SDL_BUTTON_LMASK)
          addPoint(event.motion.x, event.motion.y);
        else if(event.motion.state & SDL_BUTTON_RMASK) {
          xform = glm::rotate(basexform, (event.motion.x - initX)/400.0f, glm::vec3(0, 0, 1));
          float scale = pow(2.0, (event.motion.y - initY)/300.0f);
          // use z component to pass scale into shader so it can adjust fringe values
          xform = glm::scale(xform, glm::vec3(scale, scale, scale));
        }
      }
    }

    // draw here
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(spScreen->prog);
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(xform));

    if(points.size() > 2) {
      glVertexAttribPointer(attrPos, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), &vertex[0]); // &points[0] // &verts[0]
      glVertexAttribPointer(attrFringe, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), &fringe[0]);
      //glDrawElements(GL_TRIANGLES, triangles.size(), GL_UNSIGNED_INT, &triangles[0]);  // &elems[0]
      glDrawArrays(GL_TRIANGLES, 0, vertex.size());
    }
    checkGLError();

    // draw triangles for debugging
    if(drawTris) {
      glUseProgram(0);
      glColor4ub(255,0,0,255);
      for(unsigned int ii = 0; ii < vertex.size(); ii += 3) {
        glBegin(GL_LINE_LOOP);
        glVertex2f(vertex[ii].x, vertex[ii].y);
        glVertex2f(vertex[ii+1].x, vertex[ii+1].y);
        glVertex2f(vertex[ii+2].x, vertex[ii+2].y);
        glEnd();
      }
    }
    // end drawing

    SDL_GL_SwapWindow(window);
  }

  delete spScreen;
  //delete vbScreen;
  //delete vaScreen;

  // SDL cleanup
  SDL_GL_DeleteContext(context);
  SDL_Quit();
  return 0;
}


// stroke (variable-width)
int stroke_main(int argc, char* argv[])
{
  SDL_Init(SDL_INIT_VIDEO);

  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);  //3
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);  //2

  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow("OpenGL", 100, 100, 800, 600, SDL_WINDOW_OPENGL);
  SDL_GLContext context = SDL_GL_CreateContext(window);

  glewExperimental = GL_TRUE;
  glewInit();

  //VertexBuffer* vbScreen = new VertexBuffer(quadVertices, sizeof(quadVertices), QUAD_VERT_SIZE);
  //VertexArray* vaScreen = new VertexArray;
  std::vector<glm::vec2> points;
  std::vector<glm::vec2> normals;
  std::vector<GLfloat> miters;
  GLfloat prevx = 0, prevy = 0;

  GLfloat width = 0.001;

  std::default_random_engine rngen;
  std::uniform_int_distribution<int> rndist(1,200);
  auto rand = std::bind(rndist, rngen);

  auto addPoint = [&](int xi, int yi) {
    float x = (xi - 400.0f)/400.0f;
    float y = -(yi - 300.0f)/300.0f;
    // require minimum spacing between points
    if(std::abs(x - prevx) + std::abs(y - prevy) < 0.02)
      return;

    GLfloat w = width + rand()/10000.0f;

    int N = points.size();
    points.push_back(glm::vec2(x, y));
    points.push_back(glm::vec2(x, y));
    if(N > 0) {
      glm::vec2 n1 = glm::normalize(glm::vec2(-(y - prevy), x - prevx));
      if(N > 2) {
        normals[N-2] = glm::normalize(normals[N-2] + n1);
        normals[N-1] = normals[N-2];
      }
      else {
        normals.push_back(n1);
        normals.push_back(n1);
        miters.push_back(w);
        miters.push_back(-w);
      }
      normals.push_back(n1);
      normals.push_back(n1);
      miters.push_back(w);
      miters.push_back(-w);
    }
    prevx = x; prevy = y;
  };

  ShaderProgram* spScreen = new ShaderProgram(extrude2VS, extrude2FS);
  //vaScreen->bindAttrib(spScreen, vbScreen, "position", 2, 0);
  //vaScreen->bindAttrib(spScreen, vbScreen, "texcoord", 2, 2);

  GLint attrPos = glGetAttribLocation(spScreen->prog, "position");
  glEnableVertexAttribArray(attrPos);
  GLint attrNorm = glGetAttribLocation(spScreen->prog, "normal");
  glEnableVertexAttribArray(attrNorm);
  GLint attrMiter = glGetAttribLocation(spScreen->prog, "miter");
  glEnableVertexAttribArray(attrMiter);

  // FS uniforms
  glUniform3f(spScreen->uniform("color"), 0.0, 0.0, 1.0);
  GLuint uniModel = spScreen->uniform("model");

  //vaScreen->bind();
  spScreen->use();

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  checkGLError();

  glm::mat4 basexform;
  glm::mat4 xform;
  int initX, initY;
  SDL_Event event;
  while(1) {
    if(SDL_PollEvent(&event)) {
      if(event.type == SDL_QUIT) break;
      else if(event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_ESCAPE)
        break;
      else if(event.type == SDL_MOUSEBUTTONDOWN) {
        if(event.button.button == SDL_BUTTON_LEFT) {
          points.clear();
          normals.clear();
          miters.clear();
          prevx = 0; prevy = 0;
          xform = glm::mat4();
          basexform = glm::mat4();
          addPoint(event.button.x, event.button.y);
        }
        else if(event.button.button == SDL_BUTTON_RIGHT) {
          initX = event.button.x;
          initY = event.button.y;
          basexform = xform;
        }
      }
      else if(event.type == SDL_MOUSEMOTION) {
        if(event.motion.state & SDL_BUTTON_LMASK)
          addPoint(event.motion.x, event.motion.y);
        else if(event.motion.state & SDL_BUTTON_RMASK) {
          xform = glm::rotate(basexform, (event.motion.x - initX)/400.0f, glm::vec3(0, 0, 1));
          float scale = pow(2.0, (event.motion.y - initY)/300.0f);
          xform = glm::scale(xform, glm::vec3(scale, scale, 1.0));
        }
      }
    }

    // draw here
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(xform));

    if(points.size() > 2) {
      // attribute id, num elems, type, ?, stride (bytes), pointer
      glVertexAttribPointer(attrPos, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), &points[0]);
      glVertexAttribPointer(attrNorm, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), &normals[0]);
      glVertexAttribPointer(attrMiter, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), &miters[0]);

      //glUniform1f(uniTime, SDL_GetTicks()/1000.0);  // ms to sec
      glDrawArrays(GL_TRIANGLE_STRIP, 0, points.size());
    }
    checkGLError();
    // end drawing

    SDL_GL_SwapWindow(window);
  }

  delete spScreen;
  //delete vbScreen;
  //delete vaScreen;

  // SDL cleanup
  SDL_GL_DeleteContext(context);
  SDL_Quit();
  return 0;
}
