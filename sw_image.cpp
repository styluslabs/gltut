#include "swrender.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


namespace ImageRender {

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
Buffer2D<Color> imageTex;

class FS : public FSBase, public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) :
    gl_FragCoord(t.x + 0.5, t.y + 0.5), texCoord(baryinterp(vA.texCoord, vB.texCoord, vC.texCoord, t)) {}

  // in
  vec2 gl_FragCoord;
  vec2 texCoord;

  void main()
  {
    gl_FragColor = texture2D(imageTex, texCoord);
  }
};

} // end namespace ImageRender

class ImageSWRenderer : public SWRenderer
{
public:
  Buffer2D<Color> imageBuff;

  ImageSWRenderer(const char* imagefile)
  {
    name = "Image renderer";

    int w, h;
    //stbi_ldr_to_hdr_gamma(1.0f);  -- use w/ stbi_loadf to load as floats
    unsigned char* img = stbi_load(imagefile, &w, &h, NULL, 4);
    imageBuff = Buffer2D<Color>((Color*)img, w, h, GL_LINEAR);
  }

  void drawGlyph(LBGlyph* g) override
  {
    ImageRender::verts = &quadVertices[0];
    ImageRender::imageTex = imageBuff;

    vec2 wh = 2.0f*vec2(imageBuff.width, -imageBuff.height)/UniformBase::view_wh;
    ImageRender::bbox_center = 0.5f*wh;
    ImageRender::bbox_wh = wh;

    ShaderPipeline<ImageRender::VS, ImageRender::FS, ImageRender::Varying, BlendOver> ImagePipeline(blendOver);
    drawArrays(ImagePipeline, 6);
  }
};


#include <glm/gtc/type_ptr.hpp>

class ImageGLRenderer : public Renderer
{
public:
  const char* fillVS =
  R"(
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
  R"(
  varying vec2 texCoord;

  uniform vec4 color;
  uniform sampler2D imageTex;

  void main()
  {
    gl_FragColor = texture2D(imageTex, texCoord);
  }
  )";

  std::unique_ptr<ShaderProgram> spFill;
  GLuint imageTex;
  int imageW, imageH;
  struct { GLuint position; GLuint texcoord; } attrFill;
  struct { GLint model, scale, bbox_center, bbox_wh, color; } unifFill;

  // TODO: maybe have ShaderProgram cache the uniform, attrib indices (unordered_map?) so we don't have to
  ImageGLRenderer(const char* imagefile)
  {
    name = "GL Image";

    unsigned char* img = stbi_load(imagefile, &imageW, &imageH, NULL, 4);

    spFill.reset(new ShaderProgram(fillVS, fillFS));
    attrFill = { spFill->attrib("position_in"), spFill->attrib("texcoord_in") };
    unifFill = { spFill->uniform("model"), spFill->uniform("scale"),
        spFill->uniform("bbox_center"), spFill->uniform("bbox_wh"), spFill->uniform("color") };
    glEnableVertexAttribArray(attrFill.position);
    glEnableVertexAttribArray(attrFill.texcoord);
    glUniform1i(spFill->uniform("imageTex"), 0);

    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &imageTex);
    glBindTexture(GL_TEXTURE_2D, imageTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img);
    glBindTexture(GL_TEXTURE_2D, 0);
  }

  void drawGlyph(LBGlyph* g) override
  {
    glUseProgram(spFill->prog);
    glUniform4fv(unifFill.color, 1, glm::value_ptr(UniformBase::color));
    glUniformMatrix4fv(unifFill.model, 1, GL_FALSE, glm::value_ptr(UniformBase::model));
    glUniform1f(unifFill.scale, UniformBase::scale);
    //glUniform2f(unifFill.view_wh, UniformBase::view_wh.x, UniformBase::view_wh.y);

    vec2 wh = 2.0f*vec2(imageW, -imageH)/UniformBase::view_wh;
    glUniform2f(unifFill.bbox_center, wh.x/2, wh.y/2);
    glUniform2f(unifFill.bbox_wh, wh.x, wh.y);

    glVertexAttribPointer(attrFill.position, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][0]);
    glVertexAttribPointer(attrFill.texcoord, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GLfloat), &quadVertices[0][2]);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, imageTex);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindTexture(GL_TEXTURE_2D, 0);
  }
};
