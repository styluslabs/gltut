#include "swrender.h"
#include "stb_truetype.h"


class STBTTSWRenderer : public SWRenderer
{
public:
  stbtt_fontinfo* fontinfo;

  STBTTSWRenderer(stbtt_fontinfo* finfo) : fontinfo(finfo) { name = "STBTT"; }

  void drawGlyph(LBGlyph* g) override
  {
    // frameBuffer.data[x + frameBuffer.width*y];
    // see self->EM_per_unit in createLBFont()
    float scale = stbtt_ScaleForMappingEmToPixels(fontinfo, 400);  //stbtt_ScaleForPixelHeight(fontinfo, 400);
    float sx = UniformBase::model[0][0] * scale;
    float sy = UniformBase::model[1][1] * scale;  // sx and sy should normally be the same
    float xpos = (UniformBase::model[3][0]+1)*400;
    float ypos = (UniformBase::model[3][1]+1)*400;
    float x_shift = xpos - (float) floor(xpos);
    float y_shift = ypos - (float) floor(ypos);

    //int ascent;
    //stbtt_GetFontVMetrics(fontinfo, &ascent, 0, 0);
    //int baseline = ascent*sy;

    int x0, y0, w, h;
    unsigned char* bitmap = stbtt_GetCodepointBitmapSubpixel(fontinfo, sx, sy, x_shift, -y_shift, g->c, &w, &h, &x0, &y0);
    x0 += int(xpos);
    y0 = int(ypos) - h - y0;
    for(int ii = 0; ii < h && y0 + ii < frameBuffer.height; ++ii) {
      if(y0 + ii < 0)
        continue;
      for(int jj = 0; jj < w && x0 + jj < frameBuffer.width; ++jj) {
        if(x0 + jj < 0)
          continue;
        Color& dest = frameBuffer.data[(x0 + jj) + frameBuffer.width*(y0 + ii)];
        vec4 src = UniformBase::color;
        src.a *= bitmap[jj + w*(h - ii - 1)]/255.0f;
        for(int kk = 0; kk < 3; ++kk)
          dest[kk] = Channel(src.a*255*src[kk] + (1-src.a)*dest[kk] + 0.5f);
      }
    }
    stbtt_FreeBitmap(bitmap, NULL);
  }
};


namespace TextSum {

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
    gl_Position = model * vec4(position_in*0.5f*bbox_wh + bbox_center, 0.0f, 1.0f);
  }
};

// FS uniforms
Buffer2D<float> summedTex;

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FragCoord(t.x + 0.5, t.y + 0.5), texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)) {}

  // in
  vec2 gl_FragCoord;
  vec2 texCoord;

  float texFetchLerp(const Buffer2D<float>& texture, vec2 ij)
  {
    // texelFetch does not clamp or wrap - must do it ourselves
    //ij = clamp(ij, ivec2(0, 0), ivec2(tex_wh.x-1, tex_wh.y-1));
    float t00 = texelFetch(texture, ivec2(ij.x, ij.y), 0).r;  // implicit floor()
    float t10 = texelFetch(texture, ivec2(ij.x + 1.0f, ij.y), 0).r;
    float t01 = texelFetch(texture, ivec2(ij.x, ij.y + 1.0f), 0).r;
    float t11 = texelFetch(texture, ivec2(ij.x + 1.0f, ij.y + 1.0f), 0).r;
    vec2 f = ij - floor(ij);
    //return mix(mix(t00, t10, f.x), mix(t01, t11, f.x), f.y);
    float t0 = t00 + f.x*(t10 - t00);
    float t1 = t01 + f.x*(t11 - t01);
    return t0 + f.y*(t1 - t0);
  }

  float summedTextCov(const Buffer2D<float>& texture, vec2 st)
  {
    ivec2 tex_wh = textureSize(texture, 0);
    // for some reason, we need to shift by an extra (-0.5, -0.5) for summed case (here or in fons__getQuad)
    vec2 ij = st*vec2(tex_wh) - vec2(1.0f);
    float dx = scale*tex_wh.x/bbox_wh.x/2.0f;
    float dy = scale*tex_wh.y/bbox_wh.y/2.0f;
    float s11 = texFetchLerp(texture, ij + vec2(dx, dy));
    float s01 = texFetchLerp(texture, ij + vec2(-dx, dy));
    float s10 = texFetchLerp(texture, ij + vec2(dx,-dy));
    float s00 = texFetchLerp(texture, ij + vec2(-dx,-dy));
    float cov = (s11 - s01 - s10 + s00)/(255.0f*4.0f*dx*dy);
    return cov;
    //return (texture2D(summedTex, texCoord).r - texBiLinear(summedTex, texCoord))/100.0f + 0.5f;
    //return clamp(texFetchLerp(texture, ij)/(255.0f), 0.0, 1.0);
    //return clamp(texelFetch(texture, ivec2(ij + vec2(dx + 1.0f, dy + 1.0f)), 0).r/(48.0*48.0*255.0f), 0.0, 1.0);
  }

  void main()
  {
    //float dx = scale/bbox_wh.x/2, dy = scale/bbox_wh.y/2;
    //float s11 = texture2D(summedTex, texCoord + vec2(dx, dy)).r;
    //float s01 = texture2D(summedTex, texCoord + vec2(-dx, dy)).r;
    //float s10 = texture2D(summedTex, texCoord + vec2(dx,-dy)).r;
    //float s00 = texture2D(summedTex, texCoord + vec2(-dx,-dy)).r;
    //float cov = (s11 - s01 - s10 + s00)/(255*dx*dy*summedTex.width*summedTex.height*4);
    float cov = summedTextCov(summedTex, texCoord);

    gl_FragColor = vec4(vec3(color), cov*color.a);
  }
};

} // end namespace TextSum

class TextSumSWRenderer : public SWRenderer
{
public:
  stbtt_fontinfo* fontinfo;
  std::vector< Buffer2D<float> > atlas;
  std::vector<vec2> origins;

  TextSumSWRenderer(stbtt_fontinfo* finfo) : fontinfo(finfo)
  {
    name = "Summed coverage text";

    float s = stbtt_ScaleForMappingEmToPixels(fontinfo, 96);
    int x0, y0, w, h;
    for(int c = 32; c < 128; ++c) {
      unsigned char* bitmap = stbtt_GetCodepointBitmap(fontinfo, s, s, c, &w, &h, &x0, &y0);
      float* summed = new float[w*h];
      for(int y = 0; y < h; ++y) {
        for(int x = 0; x < w; ++x) {
          float s10 = y > 0 ? summed[(y-1)*w + x] : 0;
          float s01 = x > 0 ? summed[y*w + (x-1)] : 0;
          float s00 = x > 0 && y > 0 ? summed[(y-1)*w + (x-1)] : 0;
          summed[y*w + x] = bitmap[y*w + x] + s10 + (s01 - s00);  // /255.0f
        }
      }
      //PLATFORM_LOG("Max summed for %c: %f\n", (char)c, summed[w*h-1]);
      atlas.emplace_back(summed, w, h, GL_LINEAR); //GL_NEAREST);  //
      origins.emplace_back(x0, y0);
    }

  }

  void drawGlyph(LBGlyph* g) override
  {
    TextSum::verts = &quadVertices[0];
    TextSum::summedTex = atlas[g->c - 32];

    vec2 origin = 2.0f*(400.0f/96.0f)*origins[g->c - 32]/UniformBase::view_wh;
    vec2 wh = 2.0f*(400.0f/96.0f)*vec2(TextSum::summedTex.width, -TextSum::summedTex.height)/UniformBase::view_wh;

    TextSum::bbox_center = vec2(origin.x, -origin.y) + 0.5f*wh;
    TextSum::bbox_wh = wh;  //vec2(wh.x, -wh.y);

    ShaderPipeline<TextSum::VS, TextSum::FS, TextSum::Varying, BlendOver> TextSumPipeline(blendOver);
    drawArrays(TextSumPipeline, 6);
  }
};


// GL impl of text sum
#include <glm/gtc/type_ptr.hpp>

class TextSumGLRenderer : public Renderer
{
public:
  const char* fillVS =
  R"(#version 130
  attribute vec2 position_in;
  attribute vec2 texcoord_in;
  varying vec2 texCoord;

  uniform mat4 model;
  uniform float scale;
  uniform vec2 bbox_center;
  uniform vec2 bbox_wh;

  void main()
  {
    texCoord = texcoord_in;
    gl_Position = model * vec4(position_in*(0.5f*bbox_wh) + bbox_center, 0.0f, 1.0f);  //+ vec2(scale)
  }
  )";

  const char* fillFS =
  R"(#version 130
  varying vec2 texCoord;

  uniform vec4 color;
  uniform vec2 tex_wh;
  uniform vec2 tex_origin;
  uniform vec2 bbox_wh;
  uniform float scale;
  uniform sampler2D summedTex;

  // artifacts w/ GL_LINEAR on Intel GPU, and GLES doesn't support texture filtering for f32, so do it ourselves
  float texFetchLerp(sampler2D texture, vec2 st)
  {
    //st = tex_origin + st*tex_wh - vec2(0.5f);
    // texelFetch does not clamp or wrap - must do it ourselves
    //ist = clamp(ist, ivec2(0, 0), ivec2(tex_wh.x-1, tex_wh.y-1));
    vec2 st00 = clamp(st, tex_origin, tex_origin+tex_wh-vec2(1.0));
    vec2 st11 = clamp(st+vec2(1.0), tex_origin, tex_origin+tex_wh-vec2(1.0));
    float t00 = texelFetch(texture, ivec2(st00.x, st00.y), 0).r;  // implicit floor()
    float t10 = texelFetch(texture, ivec2(st11.x, st00.y), 0).r;
    float t01 = texelFetch(texture, ivec2(st00.x, st11.y), 0).r;
    float t11 = texelFetch(texture, ivec2(st11.x, st11.y), 0).r;
    // this seems to match behavior of texture2D() on GPU ... i.e. frac for linear interp is quantized to 8 bits
    //vec2 f = vec2(ivec2(256.0f*(st - floor(st)) + 0.5f))/256.0f;
    vec2 f = st - floor(st);
    //return mix(mix(t00, t10, f.x), mix(t01, t11, f.x), f.y);
    float t0 = t00 + f.x*(t10 - t00);
    float t1 = t01 + f.x*(t11 - t01);
    return t0 + f.y*(t1 - t0);
  }

  void main()
  {
    //ivec2 tex_wh = textureSize(summedTex, 0);  -- wrong
    // tex_wh passed in is the size of single glyph's block
    vec2 ij = tex_origin + texCoord*vec2(tex_wh) - vec2(1.0f);
    float dx = scale*tex_wh.x/bbox_wh.x/2.0f;
    float dy = scale*tex_wh.y/bbox_wh.y/2.0f;
    float s11 = texFetchLerp(summedTex, ij + vec2(dx, dy));
    float s01 = texFetchLerp(summedTex, ij + vec2(-dx, dy));
    float s10 = texFetchLerp(summedTex, ij + vec2(dx,-dy));
    float s00 = texFetchLerp(summedTex, ij + vec2(-dx,-dy));
    float cov = (s11 - s01 - s10 + s00)/(255.0*4.0*dx*dy);
    //float cov = (texture2D(summedTex, texCoord).r - texFetchLerp(summedTex, texCoord))/100.0f + 0.5f;
    gl_FragColor = vec4(vec3(color), cov*color.a);
  }
  )";

  std::unique_ptr<ShaderProgram> spFill;
  GLuint summedTex;
  struct { GLuint position; GLuint texcoord; } attrFill;
  struct { GLint model, scale, tex_wh, tex_origin, bbox_center, bbox_wh, color; } unifFill;

  stbtt_fontinfo* fontinfo;
  //std::vector< Buffer2D<float> > atlas;
  Buffer2D<float> atlas;
  std::vector<vec2> origins;

  TextSumGLRenderer(stbtt_fontinfo* finfo) : fontinfo(finfo)
  {
    name = "GL TextSum";

    int px = 40;
    float s = stbtt_ScaleForMappingEmToPixels(fontinfo, px);

    float* summed = new float[10*(px+4) * 10*(px+4)];

    int stride = 10*(px+4);
    int x0 = 2, y0 = 2;

    int xorigin, yorigin, w, h;
    for(int c = 32; c < 128; ++c) {
      unsigned char* bitmap = stbtt_GetCodepointBitmap(fontinfo, s, s, c, &w, &h, &xorigin, &yorigin);
      //float* summed = new float[w*h];
      for(int y = -2; y < px+2; ++y) {
        for(int x = -2; x < px+2; ++x) {
          float s10 = y > 0 ? summed[(y+y0-1)*stride + x+x0] : 0;
          float s01 = x > 0 ? summed[(y+y0)*stride + (x+x0-1)] : 0;
          float s00 = x > 0 && y > 0 ? summed[(y+y0-1)*stride + (x+x0-1)] : 0;
          float t11 = x >= 0 && y >= 0 && x < w && y < h ? bitmap[y*w + x] : 0;
          summed[(y+y0)*stride + x+x0] = t11 + s10 + s01 - s00;  // /255.0f
        }
      }
      x0 += px+4;
      if(x0 >= stride) {
        x0 = 2;
        y0 += px+4;
      }
      stbtt_FreeBitmap(bitmap, NULL);
      //atlas.emplace_back(summed, w, h);
      //atlas.emplace_back(bitmap, w, h);
      origins.emplace_back(xorigin, yorigin);
    }
    atlas = Buffer2D<float>(summed, stride, 10*(px+4));

    spFill.reset(new ShaderProgram(std::vector<const char*>({fillVS}), std::vector<const char*>({fillFS})));
    attrFill = { spFill->attrib("position_in"), spFill->attrib("texcoord_in") };
    unifFill = { spFill->uniform("model"), spFill->uniform("scale"), spFill->uniform("tex_wh"), spFill->uniform("tex_origin"),
        spFill->uniform("bbox_center"), spFill->uniform("bbox_wh"), spFill->uniform("color") };
    glEnableVertexAttribArray(attrFill.position);
    glEnableVertexAttribArray(attrFill.texcoord);
    glUniform1i(spFill->uniform("summedTex"), 0);

    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &summedTex);
    glBindTexture(GL_TEXTURE_2D, summedTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);  //GL_LINEAR);  //
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, atlas.width, atlas.height, 0, GL_RED, GL_FLOAT, atlas.data);

    glBindTexture(GL_TEXTURE_2D, 0);
  }

  void drawGlyph(LBGlyph* g) override
  {
    glUseProgram(spFill->prog);
    glUniform4fv(unifFill.color, 1, glm::value_ptr(UniformBase::color));
    glUniformMatrix4fv(unifFill.model, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform1f(unifFill.scale, UniformBase::scale);
    //glUniform2f(unifFill.view_wh, UniformBase::view_wh.x, UniformBase::view_wh.y);

    int cidx = g->c - 32;
    int crow = cidx / 10;
    int ccol = cidx % 10;
    int x0 = ccol * 44;
    int y0 = crow * 44;
    glUniform2f(unifFill.tex_origin, x0, y0);  //tex.width, tex.height);

    //auto& tex = atlas[g->c - 32];
    glUniform2f(unifFill.tex_wh, 44, 44);  //tex.width, tex.height);

    //vec2 wh = 20.0f*vec2(tex.width, -tex.height)/UniformBase::view_wh;
    vec2 wh = 20.0f*vec2(44, -44)/UniformBase::view_wh;
    vec2 origin = 20.0f*origins[g->c - 32]/UniformBase::view_wh;
    vec2 center = vec2(origin.x, -origin.y) + 0.5f*wh;
    glUniform2f(unifFill.bbox_center, center.x, center.y);
    glUniform2f(unifFill.bbox_wh, wh.x, wh.y);

    glVertexAttribPointer(attrFill.position, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][0]);
    glVertexAttribPointer(attrFill.texcoord, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][2]);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, summedTex);

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, tex.width, tex.height, 0, GL_RED, GL_FLOAT, tex.data);

    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
    //glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, tex.width, tex.height, 0, GL_RED, GL_UNSIGNED_BYTE, tex.data);
    //glGenerateMipmap(GL_TEXTURE_2D);

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindTexture(GL_TEXTURE_2D, 0);
  }
};
