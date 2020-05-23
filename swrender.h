#pragma once
#include "lbcommon.h"

typedef unsigned char Channel;
typedef Channel Color[4];

template<typename T>
class Buffer2D {
public:
  Buffer2D(T* d = NULL, int w = 0, int h = 0, GLenum filt = GL_NEAREST)
    : data(d), width(w), height(h), filter(filt) {}
  //Buffer2D(int w, int h) : data(new T[w*h]), width(w), height(h) {}
  //~Buffer2D() { delete data; }
  template <typename FillType = T>
  void fill(const FillType& value) { std::fill_n((FillType*)data, width*height, value); }
  template <typename SetType = T>
  void set(int x, int y, SetType val) { if(x >= 0 && x < width && y >=0 && y < height) *((SetType*)&data[x + y*width]) = val; }

  T* data;
  int width;
  int height;
  GLenum filter;
};

struct FSInput {
  FSInput() : x(0), y(0), a(0), b(0), c(0), coverage(0) {}
  FSInput& update(int _x, int _y, float _a, float _b, float _c) {
    x = _x; y = _y; a = _a; b = _b; c = _c; return *this;
  }

  // frag coords
  int x;
  int y;
  // barycentric coords
  float a;
  float b;
  float c;
  // gradients of barycentric coords
  vec2 da;
  vec2 db;
  vec2 dc;
  // experimental
  float coverage;
};

struct VaryingBase
{
  vec4 position;
};

class VSBase
{
public:
  VSBase(VaryingBase* vout, int idx) : gl_Position(vout->position) {}

  // out
  vec4& gl_Position;
};

class FSBase
{
public:
  vec4 gl_FragColor;
  //gl_FragCoord = vec4(x, y, z, 1);

  // we return a ref to allow value to be written if needed ... in OpenGL, you'd need image load/store
  template<typename T>
  T& _texelFetch(const Buffer2D<T>& texture, ivec2 ist)
  {
    // this is GL_CLAMP_TO_EDGE ... we could support other wrap configs
    ist = clamp(ist, ivec2(0, 0), ivec2(texture.width-1, texture.height-1));
    //PLATFORM_LOG("fetching (%d, %d)\n", ist.x, ist.y);
    return texture.data[ist.x + ist.y*texture.width];
  }

  vec4 texelFetch(const Buffer2D<vec4>& texture, ivec2 ist, int lod) { return _texelFetch(texture, ist); }
  vec4 texelFetch(const Buffer2D<float>& tex, ivec2 ist, int lod) { return vec4(_texelFetch(tex, ist), 0.0f, 0.0f, 1.0f); }
  vec4 texelFetch(const Buffer2D<Color>& tex, ivec2 ist, int lod)
  {
    Color& c = _texelFetch(tex, ist);
    return vec4(c[0]/255.0f, c[1]/255.0f, c[2]/255.0f, c[3]/255.0f);
  }

  template<typename T>
  vec4 texture2D(const Buffer2D<T>& texture, vec2 st)
  {
    //ivec2 ist = ivec2(floor(st*vec2(texture.width, texture.height)));
    st = st*vec2(texture.width, texture.height) - vec2(0.5f);
    vec4 t00 = texelFetch(texture, ivec2(st.x, st.y), 0);  // implicit floor()
    if(texture.filter != GL_LINEAR)
      return t00;

    vec4 t10 = texelFetch(texture, ivec2(st.x + 1, st.y), 0);
    vec4 t01 = texelFetch(texture, ivec2(st.x, st.y + 1), 0);
    vec4 t11 = texelFetch(texture, ivec2(st.x + 1, st.y + 1), 0);
    vec2 f = st - floor(st);
    //return t11*f.x*f.y + t01*(1-f.x)*f.y + t10*f.x*(1-f.y) + t00*(1-f.x)*(1-f.y);
    vec4 t0 = t00 + f.x*(t10 - t00);
    vec4 t1 = t01 + f.x*(t11 - t01);
    return t0 + f.y*(t1 - t0);
  }

  template<typename T>
  ivec2 textureSize(const Buffer2D<T>& tex, int lod) { return ivec2(tex.width, tex.height); }
};

// "typename" and "class" are 99% interchangeable in templates
template<typename T>
T baryinterp(T A, T B, T C, const FSInput& t)
{
  return t.a*A + t.b*B + t.c*C;
}

template<typename T>
T baryinterp(T A, T B, T C, float a, float b, float c)
{
  return a*A + b*B + c*C;
}

void setViewport(int left, int top, int right, int bottom);
void setDebugPixel(int x, int y);

// blending functions
class BlendOver {
public:
  Buffer2D<Color>& frameBuffer;
  BlendOver(Buffer2D<Color>& fb) : frameBuffer(fb) {}

  void blend(int x, int y, const FSBase& fs)
  {
    //vec4 src = clamp(gl_FragColor, 0.0f, 1.0f);  // expensive - remove if redundant
    const vec4& src = fs.gl_FragColor;
    // blending - see http://en.wikibooks.org/wiki/GLSL_Programming/GLUT/Transparency
    Color& dest = frameBuffer.data[x + frameBuffer.width*y];  //(fbHeight - y - 1)];
    // skip ii = 3 - destination alpha remains 1
    for(int ii = 0; ii < 3; ++ii)
      dest[ii] = Channel(src.a*255*src[ii] + (1-src.a)*dest[ii] + 0.5f);
  }
};

template<class FS>
class NoBlend {
public:
  void blend(int x, int y, const FS& fs) {}
};

// shader pipeline
class ShaderInterface
{
public:
  virtual void init(int nverts) = 0;
  virtual void run_vs(int ii) = 0;
  virtual void run_fs(int vA, int vB, int vC, const FSInput& t) const = 0;
  virtual const vec4& position(int ii) const = 0;
  virtual ~ShaderInterface() {}

  //enum ShaderType {STANDARD=0, EXACT_COVERAGE, SIGNED_DIST} shaderType;
};

template< typename VS, typename FS, typename Varying, typename Blender = NoBlend<FS> >
class ShaderPipeline : public ShaderInterface
{
public:
  NoBlend<FS> noBlend;
  Blender& blender;
  Varying* varying;

  void init(int nverts) { delete[] varying; varying = new Varying[nverts]; }
  void run_vs(int ii) { VS vs(&varying[ii], ii); vs.main(); }
  void run_fs(int vA, int vB, int vC, const FSInput& t) const
  {
    FS fs(varying[vA], varying[vB], varying[vC], t);
    fs.main();
    blender.blend(t.x, t.y, fs);
  }
  const vec4& position(int ii) const { return varying[ii].position; }

  ShaderPipeline(Blender& b) : blender(b), varying(NULL) {}
  ShaderPipeline() : blender(noBlend), varying(NULL) {}
  ~ShaderPipeline() { delete[] varying; }
};

void drawArrays(ShaderInterface& prog, int nverts, int geommode = GL_TRIANGLES);
void drawElements(ShaderInterface& prog, int* indices, int nindices, int nverts, int geommode = GL_TRIANGLES);

// abstract Renderer class for rendering buffer (texture) to screen
#include <memory>

class SWRenderer : public Renderer
{
public:
  Buffer2D<Color> frameBuffer;  // we could make this static if we want to save memory
  BlendOver blendOver;

  std::unique_ptr<ShaderProgram> spScreen;
  //GLuint attrPos_Screen, attrTex_Screen;
  //GLuint vao;
  GLuint texScreen;

  static unsigned int clearColor;

  SWRenderer();
  ~SWRenderer() override;
  void beginFrame() override;
  void endFrame() override;
  bool isSWRenderer() const override { return true; }
};
