#include <string>
#include <chrono>

#define GLPP_IMPLEMENTATION
#include "glpp.h"

#include "swrender.h"
#include "lbfont.h"

#include <glm/gtc/matrix_transform.hpp>

// SDL #defines main in Windows!
#include <SDL.h>
#ifdef main
#undef main
#else
#define SDL_main main
#endif

//#include "gl_lb.cpp"
#include "gl_winding.cpp"
#include "sw_lb.cpp"
//#include "sw_exact.cpp"
#include "sw_sdf.cpp"
#include "sw_winding.cpp"
#include "sw_stbtt.cpp"
#include "sw_stroke.cpp"
#include "sw_coverage.cpp"
#include "sw_trapezoid.cpp"
#include "sw_image.cpp"

static double calcBrightness(Buffer2D<Color>& buff)
{
  double brightness = 0;
  for(int y = 0; y < buff.height; ++y) {
    for(int x = 0; x < buff.width; ++x) {
      //0.2126*red() + 0.7152*green() + 0.0722*blue()
      // just use green component, assuming greyscale only; we apply alpha assuming white BG
      Color& src = buff.data[y*buff.width + x];
      double a = src[3]/255.0;
      brightness += 255 - (a*src[1] + (1 - a)*255 + 0.5); //255 - buff.data[y*buff.width + x][1];
    }
  }
  return brightness;
}

// render with OpenGL
int SDL_main(int argc, char* argv[])
{
  SDL_Init(SDL_INIT_VIDEO);

  ShaderProgram::headerVS = "#version 120\n";
  R"(#version 300 es
  #define attribute in
  #define varying out
  )";

  ShaderProgram::headerFS = "#version 120\n";
  R"(#version 300 es
  precision highp float;
  #define varying in

  #define gl_FragColor fragColor
  out vec4 fragColor;
  )";

  // profiles (core vs. compat) are only meaningful for 3.2+; core requires using VBOs and VAOs
  // vmware linux supports 3.3 core and compat (and ES3.0); vmware windows only supports 3.0 and 3.3 core
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);  //SDL_GL_CONTEXT_PROFILE_ES
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);

  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

  SDL_Window* window = SDL_CreateWindow("LBFont - GL Renderer", 100, 100, 800, 800, SDL_WINDOW_OPENGL);
  SDL_GLContext context = SDL_GL_CreateContext(window);

  glewExperimental = GL_TRUE;
  glewInit();

  const char* fontfile = argc > 2 ? argv[2] : "../../nanovg/example/Roboto-Regular.ttf";
  LBFont* font = createLBFont(fontfile);
  if(!font) {
    PLATFORM_LOG("Error loading font from: %s\n", fontfile);
    return -1;
  }

  // fill color
  vec4 color(0.0f, 0.0f, 0.0f, 0.5f);

  // glyph for scribble shape
  LBGlyph* scribble = NULL;
  std::vector<stbtt_vertex> scribble_verts;
  size_t subpath_start = 0;

  // TODO: use uniform buffer object to share uniforms between shaders!
  { // scope so renderers are deleted before GL context
  //LBGLRenderer lbGLRenderer;
  WindingGLRenderer windingGLRenderer;
  LBSWRenderer lbSWRenderer;
  //ExactSWRenderer exactSWRenderer;
  SDFSWRenderer sdfSWRenderer;
  WindingSWRenderer windingSWRenderer;
  StrokeSWRenderer strokeSWRenderer;
  STBTTSWRenderer stbttSWRenderer(font->fontinfo);
  CoverageSWRenderer coverageSWRenderer;
  TrapezoidSWRenderer trapezoidSWRenderer;
  //StrokeSWRenderer strokeSWRenderer;
  TextSumSWRenderer textSumSWRenderer(font->fontinfo);
  TextSumGLRenderer textSumGLRenderer(font->fontinfo);
  ImageSWRenderer imageSWRenderer("/home/mwhite/graphics/gltut/sample.png");
  ImageGLRenderer imageGLRenderer("/home/mwhite/graphics/gltut/sample.png");

  glEnable(GL_BLEND);
  GL_CHECK(glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));
  // opengl timer ... doesn't seem to return correct values in Linux VM
  //GLuint timerQuery;
  //glGenQueries(1, &timerQuery);
  //uint64_t timerAccum = 0;

  setViewport(0, 0, 800, 800);  // SW rasterizer setup
  float aastep_fill = 0.006f;

  Renderer* currRenderer = &windingGLRenderer;
  mat4 basexform(1.0f), xform(1.0f);
  int initX = 0, initY = 0;
  int mouseX = 0, mouseY = 0;
  SDL_Keymod keymod = KMOD_NONE;
  SDL_Event event;
  //bool drawTris = false;
  std::string teststring("WXYZIALMTV");
  // for FPS counter
  auto start = std::chrono::high_resolution_clock::now();
  int nframes = 0;
  bool run = true;
  bool dirty = true;
  int nlayers = 0;
  bool analyze = false;
  bool contFPS = false;
  bool drawmode = false;
  bool whiteOnBlack = false;
  double prevbrightness = 0;
  // event loop
  while(run) {
    // use while() to empty event queue before rendering frame, so that user's most recent input is accounted for
    while((dirty || contFPS) ? SDL_PollEvent(&event) : SDL_WaitEvent(&event)) {
      if(event.type == SDL_QUIT)
        run = false;
      else if(event.type == SDL_KEYDOWN) {
        if(event.key.keysym.sym == SDLK_ESCAPE)
          run = false;
        // renderer seletion
        else if(event.key.keysym.sym == SDLK_F1)
          currRenderer = &windingGLRenderer;
        else if(event.key.keysym.sym == SDLK_F2)
          currRenderer = &lbSWRenderer;
        else if(event.key.keysym.sym == SDLK_F3)
          currRenderer = &textSumGLRenderer;  //&imageGLRenderer;  //&sdfSWRenderer;
        else if(event.key.keysym.sym == SDLK_F4)
          currRenderer = &textSumSWRenderer;  //&imageSWRenderer;  //
        else if(event.key.keysym.sym == SDLK_F5)
          currRenderer = &windingSWRenderer;  //&exactSWRenderer;
        else if(event.key.keysym.sym == SDLK_F6)
          currRenderer = &trapezoidSWRenderer;
        else if(event.key.keysym.sym == SDLK_F7)
          currRenderer = &stbttSWRenderer;
        else if(event.key.keysym.sym == SDLK_F8)
          currRenderer = &coverageSWRenderer;
        // misc config
        else if(event.key.keysym.sym == SDLK_F9) {
          whiteOnBlack = !whiteOnBlack;
          color = whiteOnBlack ? vec4(1.0f, 1.0f, 1.0f, 1.0f) : vec4(0.0f, 0.0f, 0.0f, 1.0f);
        }
          //imageSWRenderer.imageBuff.filter = imageSWRenderer.imageBuff.filter == GL_LINEAR ? GL_NEAREST : GL_LINEAR;
          //SDFGen::drawTris = !SDFGen::drawTris;
        else if(event.key.keysym.sym == SDLK_F10) {
          nlayers = (event.key.keysym.mod & KMOD_SHIFT) ? 2*nlayers : !nlayers;
          PLATFORM_LOG("nlayers: %d\n", nlayers);
        }
        else if(event.key.keysym.sym == SDLK_F11)
          analyze = !analyze;
        else if(event.key.keysym.sym == SDLK_F12)
          contFPS = !contFPS;
        else if(event.key.keysym.sym == SDLK_INSERT) {
          xform = mat4(1.0f);
          drawmode = !drawmode;
        }
        else if(event.key.keysym.sym == SDLK_BACKSPACE)
          teststring.pop_back();
        else if(event.key.keysym.sym == SDLK_DELETE) {
          teststring.clear();
          scribble_verts.clear();
          subpath_start = 0;
          freeLBGlyph(scribble);
          scribble = NULL;
        }
        else if(event.key.keysym.sym == SDLK_HOME) {
          // reset transform
          basexform = mat4(1.0f);
          xform = mat4(1.0f);
        }
        else if(event.key.keysym.sym == SDLK_EQUALS || event.key.keysym.sym == SDLK_MINUS) {
          int dir = event.key.keysym.sym == SDLK_EQUALS ? 1 : -1;
          if(event.key.keysym.mod & KMOD_ALT)
            color.a += 0.1f*dir;
          else if(event.key.keysym.mod & KMOD_CTRL)
            xform = glm::scale(mat4(1.0f), vec3(dir > 0 ? 1.25f : 0.8f))*xform;
        }
      }
      else if(event.type == SDL_TEXTINPUT) {
        teststring.append(event.text.text);

        char c = teststring.back();
        if(c <= font->num_glyphs) {
          LBGlyph* g = font->glyphs[int(c)];
          float bboxarea = (g->x_max - g->x_min)*(g->y_max - g->y_min);
          PLATFORM_LOG("%c : bbox/area = %f/%f = %f\n", c, bboxarea, g->area, bboxarea/g->area);
        }
      }
      else if(event.type == SDL_MOUSEBUTTONDOWN) {
        if(drawmode) {
          unsigned char cmd = (scribble_verts.empty() || event.button.button == SDL_BUTTON_RIGHT) ? STBTT_vmove : STBTT_vline;
          if(cmd == STBTT_vmove && !scribble_verts.empty()) {
            scribble_verts.push_back(scribble_verts[subpath_start]);
            scribble_verts.back().type = STBTT_vline;
            subpath_start = scribble_verts.size();
          }
          stbtt_vertex v = {short(event.button.x - 400), short(400 - event.button.y), 0, 0, 0, 0, cmd, 0};
          scribble_verts.push_back(v);
          if(scribble_verts.size() > subpath_start + 2) {
            freeLBGlyph(scribble);
            scribble = (LBGlyph*)calloc(1, sizeof(LBGlyph));
            scribble_verts.push_back(scribble_verts[subpath_start]);
            scribble_verts.back().type = STBTT_vline;
            initLBGlyphPts(scribble, 1/400.0f, scribble_verts.size(), scribble_verts.data());
            scribble_verts.pop_back();
          }
        }
        else if(event.button.button == SDL_BUTTON_LEFT) {
          keymod = SDL_GetModState();
          initX = event.button.x;
          initY = event.button.y;
          basexform = xform;
        }
        else if(event.button.button == SDL_BUTTON_RIGHT) {
          // we'll probably find a better use for right button, in which case add a modifier for debug select
          PLATFORM_LOG("Setting debug pixel (%d,%d)\n", event.button.x, 800 - event.button.y);
          setDebugPixel(event.button.x, event.button.y);
        }
      }
      else if(event.type == SDL_MOUSEMOTION) {
        if(drawmode)
          continue;
        mouseX = event.motion.x;  mouseY = event.motion.y;
        float fine = (keymod & KMOD_ALT) ? 0.1f : 1.0f;
        if(event.motion.state & SDL_BUTTON_LMASK) {
          if(keymod & KMOD_SHIFT) {
            xform = glm::rotate(mat4(1.0f), fine*(event.motion.x - initX)/400.0f, vec3(0, 0, 1))*basexform;
          }
          else if(keymod & KMOD_CTRL) {
            float scale = pow(2.0f, fine*(event.motion.x - initX)/400.0f);
            // use z component to pass scale into shader so it can adjust fringe values
            xform = glm::scale(mat4(1.0f), vec3(scale))*basexform;
          }
          else {
            xform = glm::translate(mat4(1.0f),
                vec3(fine*(event.motion.x - initX)/400.0f, -fine*(event.motion.y - initY)/400.0f, 0))*basexform;
          }
        }
      }
      else if(event.type == SDL_WINDOWEVENT) {
        if(event.window.event != SDL_WINDOWEVENT_EXPOSED && event.window.event != SDL_WINDOWEVENT_MOVED
            && event.window.event != SDL_WINDOWEVENT_RESIZED)
          continue;
      }
      else
        continue;
      dirty = true;
    }
    dirty = false;
    std::string title = std::string(currRenderer->name) + (drawmode ? " - Draw mode" : "");
    if(title != SDL_GetWindowTitle(window))
      SDL_SetWindowTitle(window, title.c_str());

    // set common uniforms
    UniformBase::model = xform;
    UniformBase::scale = 1/(400.0f*xform[2][2]);
    UniformBase::view_wh = vec2(800, 800);
    UniformBase::aastep = aastep_fill;
    UniformBase::color = color;

    // begin drawing
    SWRenderer::clearColor = whiteOnBlack ? 0xFF000000 : 0xFFFFFFFF;
    whiteOnBlack ? glClearColor(0.0f, 0.0f, 0.0f, 1.0f) : glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glBeginQuery(GL_TIME_ELAPSED, timerQuery);
    currRenderer->beginFrame();

    // render scribble path if not empty
    if(scribble && scribble->num_verts > 0)
      currRenderer->drawGlyph(scribble);

    // option to write multiple strings to test overlapping paths
    std::vector<std::string> ss = {teststring};
    for(int ii = 0; ii < nlayers; ++ii)
      ss.emplace_back(ss.back().rbegin(), ss.back().rend());
    for(const std::string& s : ss) {
      mat4 textxform = xform;
      for(size_t ii = 0; ii < s.size(); ++ii) {
        char c = s[ii];
        if(c > font->num_glyphs) continue;
        LBGlyph* g = font->glyphs[int(c)];
        if(!g) continue;
        UniformBase::model = textxform;
        // whitespace characters will have no vertices
        if(g->num_verts > 0)
          currRenderer->drawGlyph(g);
        textxform = glm::translate(textxform, vec3(g->advance, 0, 0));
      }
    }

    if(currRenderer->isSWRenderer()) {
      Buffer2D<Color>& fb = static_cast<SWRenderer*>(currRenderer)->frameBuffer;
      if(analyze) {
        double brightness = calcBrightness(fb);
        PLATFORM_LOG("Total brightness: %f; delta: %f\n", brightness, brightness - prevbrightness);
        prevbrightness = brightness;
      }
      // draw mouse position as a red pixel so it can seen in magnifer
      fb.set(mouseX, fb.height - mouseY, 0xFF0000FF);
    }
    currRenderer->endFrame();
    //glEndQuery(GL_TIME_ELAPSED);
    GL_CHECK("End of frame");
    SDL_GL_SwapWindow(window);
    // end drawing
    // note that GL_QUERY_RESULT blocks until result is available, so if we want to use for a subset of calls
    //  we should defer until endFrame() (or just before glBeginQuery of next frame)
    // see http://ephenationopengl.blogspot.com/2012/01/measuring-graphics-performance.html
    //GLuint timerQueryResult;
    //glGetQueryObjectuiv(timerQuery, GL_QUERY_RESULT, &timerQueryResult);
    //timerAccum += timerQueryResult;

    nframes++;
    auto now = std::chrono::high_resolution_clock::now();
    float seconds = std::chrono::duration<float>(now - start).count();
    if(contFPS && nframes > 1 && seconds > 1) {
      //PLATFORM_LOG("FPS (%s): %f ; GL FPS: %f\n", rendererName, nframes/seconds, nframes/(timerAccum/1E9));
      PLATFORM_LOG("FPS (%s): %f\n", currRenderer->name, nframes/seconds);
      start = now;
      nframes = 0;
      //timerAccum = 0;
    }
  }

  //glDeleteQueries(1, &timerQuery);
  } // delete renderers
  freeLBGlyph(scribble);
  freeLBFont(font);

  // SDL cleanup
  SDL_GL_DeleteContext(context);
  SDL_Quit();
  return 0;
}
