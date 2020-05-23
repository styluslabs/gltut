#pragma once
#include <iostream>
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>

#ifdef __ANDROID__
#include <android/log.h>
#define PLATFORM_LOG(...) __android_log_print(ANDROID_LOG_VERBOSE, "GL++",  __VA_ARGS__);
#elif defined(_WIN32)
// only downside to SDL_Log is that messages are limited to 4KB
#include <SDL_log.h>
#define PLATFORM_LOG SDL_Log
#else
// note that fprint(stdout, ...) requires fflush for Qt Creator (but not stderr)
#define PLATFORM_LOG(...) fprintf(stderr, __VA_ARGS__)
#endif

#ifndef NDEBUG
void _checkGLError(const char* call, const char *file, int line);
#define checkGLError(s) _checkGLError(#s, __FILE__, __LINE__)
#define GL_CHECK(x) x; checkGLError(#x)
#else
#define GL_CHECK(x) x
#endif

class ShaderProgram {
public:
  bool ok;
  GLuint prog;
  GLuint vert;
  GLuint frag;

  ShaderProgram(const std::vector<const char*>& s_vert, const std::vector<const char*>& s_frag) : ok(true)
  {
    loadProgram(s_vert, s_frag);
    GL_CHECK(glUseProgram(prog));
    currProg = prog;
  }

  ShaderProgram(const char* s_vert, const char* s_frag) : ShaderProgram({headerVS, s_vert}, {headerFS, s_frag,}) {}

  ~ShaderProgram()
  {
    GL_CHECK(glDeleteProgram(prog));
    GL_CHECK(glDeleteShader(frag));
    GL_CHECK(glDeleteShader(vert));
  }

  void loadProgram(const std::vector<const char*>& s_vert, const std::vector<const char*>& s_frag)
  {
    vert = loadShader(s_vert, GL_VERTEX_SHADER);
    frag = loadShader(s_frag, GL_FRAGMENT_SHADER);
    prog = GL_CHECK(glCreateProgram());
    GL_CHECK(glAttachShader(prog, vert));
    GL_CHECK(glAttachShader(prog, frag));
    GL_CHECK(glLinkProgram(prog));
  }

  void use()
  {
    if(prog != currProg) {
      GL_CHECK(glUseProgram(prog));
      currProg = prog;
    }
  }

  static void useProgram(GLuint p)
  {
    GL_CHECK(glUseProgram(p));
    currProg = p;
  }

  // previously tried caching results in unordered_map, but lookup turned out to be slow
  GLint uniform(const char* s) { return glGetUniformLocation(prog, s); }
  GLuint attrib(const char* s) { return glGetAttribLocation(prog, s); }

  static const char* headerVS;
  static const char* headerFS;

private:
  GLuint loadShader(const std::vector<const char*>& src, GLuint type)
  {
    GLint status;
    GLuint id = glCreateShader(type);
    glShaderSource(id, src.size(), src.data(), NULL);
    //glShaderSource(id, 1, &src, NULL);
    glCompileShader(id);
    glGetShaderiv(id, GL_COMPILE_STATUS, &status);
    if(status != GL_TRUE) {
      char buffer[512];
      glGetShaderInfoLog(id, 512, NULL, buffer);
      PLATFORM_LOG("Error compiling shader:\n");
      for(const char* s : src)
        PLATFORM_LOG(s);
      PLATFORM_LOG(buffer);
      ok = false;
    }
    return id;
  }

  static GLuint currProg;
};

class VertexBuffer
{
public:
  GLuint VBO;
  size_t stride;

  VertexBuffer(void* buffer, size_t size, size_t _stride) : stride(_stride)
  {
    GL_CHECK(glGenBuffers(1, &VBO));
    GL_CHECK(glBindBuffer(GL_ARRAY_BUFFER, VBO));
    GL_CHECK(glBufferData(GL_ARRAY_BUFFER, size, buffer, GL_STATIC_DRAW));
    currVBO = VBO;
  }

  ~VertexBuffer()
  {
    GL_CHECK(glDeleteBuffers(1, &VBO));
  }

  void bind() const
  {
    if(VBO != currVBO) {
      GL_CHECK(glBindBuffer(GL_ARRAY_BUFFER, VBO));
      currVBO = VBO;
    }
  }

private:
  static GLuint currVBO;
};

class VertexArray
{
public:
  GLuint VAO;

  VertexArray()
  {
    GL_CHECK(glGenVertexArrays(1, &VAO));
    GL_CHECK(glBindVertexArray(VAO));
    currVAO = VAO;
  }

  ~VertexArray()
  {
    GL_CHECK(glDeleteVertexArrays(1, &VAO));
  }

  void bind() const
  {
    if(VAO != currVAO) {
      GL_CHECK(glBindVertexArray(VAO));
      currVAO = VAO;
    }
  }

  void bindAttrib(const ShaderProgram* prog, const VertexBuffer* buff, const char* attrib, GLuint nelem, GLuint offset)
  {
    bind();
    buff->bind();

    GLint id = GL_CHECK(glGetAttribLocation(prog->prog, attrib));
    GL_CHECK(glEnableVertexAttribArray(id));
    // glVertexAttribPointer(attribute, num elements (1-4), element type, normalized, stride, offset)
    GL_CHECK(glVertexAttribPointer(id, nelem, GL_FLOAT, GL_FALSE, (buff->stride*sizeof(GLfloat)), (void*)(offset*sizeof(GLfloat))));
  }

private:
  static GLuint currVAO;
};

//#define QUAD_VERT_SIZE 4
// {x, y, u, v}, where u,v are texture coords
extern GLfloat quadVertices[6][4];

#ifdef GLPP_IMPLEMENTATION

GLuint VertexBuffer::currVBO = -1;
GLuint VertexArray::currVAO = -1;
GLuint ShaderProgram::currProg = -1;
const char* ShaderProgram::headerVS = "";
const char* ShaderProgram::headerFS = "";

GLfloat quadVertices[][4] = {
  {-1.0f,  1.0f,  0.0f, 1.0f},
  { 1.0f,  1.0f,  1.0f, 1.0f},
  { 1.0f, -1.0f,  1.0f, 0.0f},

  { 1.0f, -1.0f,  1.0f, 0.0f},
  {-1.0f, -1.0f,  0.0f, 0.0f},
  {-1.0f,  1.0f,  0.0f, 1.0f}
};

#ifndef NDEBUG
void _checkGLError(const char* call, const char* file, int line) {
  GLenum err(glGetError());
  while(err != GL_NO_ERROR) {
    const char* error;
    switch(err) {
      case GL_INVALID_OPERATION:      error="INVALID_OPERATION";      break;
      case GL_INVALID_ENUM:           error="INVALID_ENUM";           break;
      case GL_INVALID_VALUE:          error="INVALID_VALUE";          break;
      case GL_OUT_OF_MEMORY:          error="OUT_OF_MEMORY";          break;
      case GL_INVALID_FRAMEBUFFER_OPERATION:  error="INVALID_FRAMEBUFFER_OPERATION";  break;
    }
    PLATFORM_LOG("GL_%s - %s at %s:%d\n", error, call, file, line);
    err = glGetError();
  }
}
#endif

#endif
