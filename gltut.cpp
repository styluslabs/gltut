#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

// Eigen is a popular alternative to GLM for more advanced linear algebra
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "glpp.h"

#define CUBE_VERT_SIZE 8
GLfloat cubeVertices[] = {
  // X      Y     Z     R     G     B     U     V
  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
   0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
   0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
   0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
  -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

  -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
   0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
  -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
  -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

  -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
  -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
  -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
  -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,

   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
   0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
   0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
   0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
   0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,

  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
   0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
   0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
   0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
  -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
  -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,

  -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,
   0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
   0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f,
  -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,
  -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f,

  -1.0f, -1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
   1.0f, -1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
   1.0f,  1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f,
   1.0f,  1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f,
  -1.0f,  1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
  -1.0f, -1.0f, -0.5f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
};

#define QUAD_VERT_SIZE 4
GLfloat quadVertices[] = {
  -1.0f,  1.0f,  0.0f, 1.0f,
   1.0f,  1.0f,  1.0f, 1.0f,
   1.0f, -1.0f,  1.0f, 0.0f,

   1.0f, -1.0f,  1.0f, 0.0f,
  -1.0f, -1.0f,  0.0f, 0.0f,
  -1.0f,  1.0f,  0.0f, 1.0f
};

// define vertices constituting our triangles
GLuint elements[] = {
  0, 1, 2,
  2, 3, 0
};

const char* sceneVert =
R"(#version 120

attribute vec3 position;
attribute vec3 color;
attribute vec2 texcoord;

varying vec3 Color;
varying vec2 Texcoord;

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;
uniform vec3 overrideColor;

void main()
{
  Color = overrideColor * color;
  Texcoord = texcoord;
  gl_Position = proj * view * model * vec4(position, 1.0);
}
)";

const char* sceneFrag =
R"(#version 120

varying vec3 Color;
varying vec2 Texcoord;

uniform sampler2D texKitten;
uniform sampler2D texPuppy;

void main()
{
  vec4 colKitten = texture2D(texKitten, Texcoord);
  vec4 colPuppy = texture2D(texPuppy, Texcoord);
  gl_FragColor = mix(colKitten, colPuppy, 0.5) * vec4(Color, 1.0);
}
)";

const char* screenVert =
R"(#version 120

attribute vec2 position;
attribute vec2 texcoord;

varying vec2 Texcoord;

void main()
{
   Texcoord = texcoord;
   gl_Position = vec4(position, 0.0, 1.0);
}
)";

const char* screenFrag =
R"(#version 120

varying vec2 Texcoord;
uniform sampler2D texFramebuffer;

const float blurSizeH = 1.0 / 300.0;
const float blurSizeV = 1.0 / 200.0;
void main()
{
  vec4 sum = vec4(0.0);
  for (int x = -4; x <= 4; x++) {
    for (int y = -4; y <= 4; y++)
      sum += texture2D(texFramebuffer, vec2(Texcoord.x + x * blurSizeH, Texcoord.y + y * blurSizeV)) / 81.0;
  }
  gl_FragColor = sum;
  // gl_FragColor = texture2D(texFramebuffer, Texcoord);
}
)";

GLuint loadTexture(const char* path)
{
  GLuint texbuf;
  GL_CHECK(glGenTextures(1, &texbuf));
  GL_CHECK(glBindTexture(GL_TEXTURE_2D, texbuf));

  SDL_Surface* teximage = IMG_Load(path);
  if(!teximage) {
    std::cout << "error loading image!";
    return texbuf;
  }
  GL_CHECK(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, teximage->w, teximage->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, teximage->pixels));
  SDL_FreeSurface(teximage);

  GL_CHECK(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE));
  GL_CHECK(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE));
  GL_CHECK(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR));
  GL_CHECK(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));
  return texbuf;
}

// saveFramebuff(800, 600, "/home/mwhite/maps/gltut/img.bmp");
void saveFramebuff(int w, int h, const char* path)
{
  GLubyte* pixels = new GLubyte[w*h*4];
  glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
  Uint32 rmask = 0xff000000, gmask = 0x00ff0000, bmask = 0x0000ff00, amask = 0x000000ff;
#else
  Uint32 rmask = 0x000000ff, gmask = 0x0000ff00, bmask = 0x00ff0000, amask = 0xff000000;
#endif
  SDL_Surface* saveSurface = SDL_CreateRGBSurfaceFrom(pixels, w, h, 32, w * 4, rmask, gmask, bmask, amask);
  SDL_SaveBMP(saveSurface, path);
  delete[] pixels;
}

// Ref: tutorial from https://open.gl/
int main(int argc, char* argv[])
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

  VertexBuffer* vbScene = new VertexBuffer(cubeVertices, sizeof(cubeVertices), CUBE_VERT_SIZE);
  VertexArray* vaScene = new VertexArray;

  VertexBuffer* vbScreen = new VertexBuffer(quadVertices, sizeof(quadVertices), QUAD_VERT_SIZE);
  VertexArray* vaScreen = new VertexArray;

  // element buffer
  /*GLuint ebo;
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
  */

  // Load textures
  GLuint texKitten = loadTexture("/home/mwhite/graphics/gltut/sample.png");
  GLuint texPuppy = loadTexture("/home/mwhite/graphics/gltut/sample2.png");

  // load shaders
  ShaderProgram* spScene = new ShaderProgram(sceneVert, sceneFrag);
  // bind shader attributes to vertex array elements
  vaScene->bindAttrib(spScene, vbScene, "position", 3, 0);
  vaScene->bindAttrib(spScene, vbScene, "color", 3, 3);
  vaScene->bindAttrib(spScene, vbScene, "texcoord", 2, 6);

  glUniform1i(spScene->uniform("texKitten"), 0);  // GL_TEXTURE0
  glUniform1i(spScene->uniform("texPuppy"), 1);  // GL_TEXTURE1

  GLint uniModel = spScene->uniform("model");
  GLint uniColor = spScene->uniform("overrideColor");

  // Set up projection
  glm::mat4 view = glm::lookAt(
    glm::vec3(2.5f, 2.5f, 2.0f),
    glm::vec3(0.0f, 0.0f, 0.0f),
    glm::vec3(0.0f, 0.0f, 1.0f)
  );
  glUniformMatrix4fv(spScene->uniform("view"), 1, GL_FALSE, glm::value_ptr(view));

  glm::mat4 proj = glm::perspective(45.0f, 800.0f / 600.0f, 1.0f, 10.0f);
  glUniformMatrix4fv(spScene->uniform("proj"), 1, GL_FALSE, glm::value_ptr(proj));

  // screen shader
  ShaderProgram* spScreen = new ShaderProgram(screenVert, screenFrag);
  vaScreen->bindAttrib(spScreen, vbScreen, "position", 2, 0);
  vaScreen->bindAttrib(spScreen, vbScreen, "texcoord", 2, 2);

  glUniform1i(spScreen->uniform("texFramebuffer"), 0);

  // create framebuffer into which we will render scene
  GLuint frameBuffer;
  glGenFramebuffers(1, &frameBuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

  // Create texture to hold color buffer
  GLuint texColorBuffer;
  glGenTextures(1, &texColorBuffer);
  glBindTexture(GL_TEXTURE_2D, texColorBuffer);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);

  // Create Renderbuffer Object to hold depth and stencil buffers
  GLuint rboDepthStencil;
  glGenRenderbuffers(1, &rboDepthStencil);
  glBindRenderbuffer(GL_RENDERBUFFER, rboDepthStencil);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 800, 600);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rboDepthStencil);

  if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
  {
    std::cout << "Error creating framebuffer\n";
  }
  checkGLError();

  SDL_Event windowEvent;
  while(1) {
    if(SDL_PollEvent(&windowEvent)) {
      if (windowEvent.type == SDL_QUIT) break;
      else if(windowEvent.type == SDL_KEYUP && windowEvent.key.keysym.sym == SDLK_ESCAPE) break;
    }

    // draw here

    // Bind our framebuffer and draw 3D scene (spinning cube)
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    vaScene->bind();
    glEnable(GL_DEPTH_TEST);
    spScene->use();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texKitten);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texPuppy);

    // Clear the screen to white
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // model transformation
    glm::mat4 trans;
    trans = glm::rotate(trans, (float)clock() / (float)CLOCKS_PER_SEC * 5.0f, glm::vec3(0.0f, 0.0f, 1.0f));
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(trans));

    //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glDrawArrays(GL_TRIANGLES, 0, 36);

    // how stencil buffer is used here: we write 0xFF when drawing floor, then when we draw cube, we mask (reject)
    //  any pixels which don't have 0xFF in stencil buffer - this limits reflection to the floor
    glEnable(GL_STENCIL_TEST);
    // Draw floor
    glStencilFunc(GL_ALWAYS, 1, 0xFF);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    glStencilMask(0xFF);
    glDepthMask(GL_FALSE);  // don't include floor in depth buffer so that we can draw reflected cube
    glClear(GL_STENCIL_BUFFER_BIT);

    glDrawArrays(GL_TRIANGLES, 36, 6);

    // draw reflected, dimmed version of cube to simulate reflection
    glStencilFunc(GL_EQUAL, 1, 0xFF);
    glStencilMask(0x00);
    glDepthMask(GL_TRUE);

    trans = glm::scale(glm::translate(trans, glm::vec3(0, 0, -1)), glm::vec3(1, 1, -1));
    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(trans));

    glUniform3f(uniColor, 0.3f, 0.3f, 0.3f);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glUniform3f(uniColor, 1.0f, 1.0f, 1.0f);
    glDisable(GL_STENCIL_TEST);

    // Bind default framebuffer and draw contents of our framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    vaScreen->bind();
    glDisable(GL_DEPTH_TEST);
    spScreen->use();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texColorBuffer);

    glDrawArrays(GL_TRIANGLES, 0, 6);

    // end drawing

    SDL_GL_SwapWindow(window);
  }

  glDeleteRenderbuffers(1, &rboDepthStencil);
  glDeleteTextures(1, &texColorBuffer);
  glDeleteFramebuffers(1, &frameBuffer);

  // GL cleanup
  glDeleteTextures(1, &texKitten);
  glDeleteTextures(1, &texPuppy);

  delete spScene;
  delete spScreen;
  delete vbScene;
  delete vaScene;
  delete vbScreen;
  delete vaScreen;

  // SDL cleanup
  SDL_GL_DeleteContext(context);
  SDL_Quit();
  return 0;
}
