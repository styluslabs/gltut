#include "swrender.h"

// shader for exact coverage
// note: namespace must be named - anonymous namespace will allow enclosed "using namespace ..." to break out
namespace ExactFill {

using namespace glm;

struct Varying : public VaryingBase {};

const Vertex* verts;

class VS : public VSBase, public UniformBase {
public:
  VS(Varying* vout, int idx) : VSBase(vout, idx), position((float)verts[idx].x, (float)verts[idx].y) {}

  // in
  vec2 position;

  void main()
  {
    gl_Position = model * vec4(position, 0.0f, 1.0f);
  }
};

Buffer2D<Color> frameBuffer;

class FS : public UniformBase {
public:
  FS(const Varying& vA, const Varying& vB, const Varying& vC, const FSInput& t) { run(t.x, t.y, t.coverage); }

  void main() {}
  void run(int x, int y, float coverage)
  {
    Color& dest = frameBuffer.data[x + frameBuffer.width*y];
    dest[3] = (unsigned char)min(color.a*255.0f, dest[3] + coverage*color.a*255.0f + 0.5f);
    for(int ii = 0; ii < 3; ++ii)
      dest[ii] = (unsigned char)(color[ii]);
  }
};

} // end namespace ExactFill


class ExactSWRenderer : public SWRenderer
{
public:
  ExactSWRenderer() { ExactFill::frameBuffer = frameBuffer;  name = "SW Exact"; }

  void beginFrame() override
  {
    frameBuffer.fill<unsigned int>(0x00FFFFFF);
  }

  void drawGlyph(LBGlyph* g) override
  {
    // vertices
    ExactFill::verts = g->triangle_vertices;

    ShaderPipeline<ExactFill::VS, ExactFill::FS, ExactFill::Varying> ExactFillPipeline(ShaderInterface::EXACT_COVERAGE);
    drawArrays(ExactFillPipeline, g->num_triangle_vertices);
  }
};
