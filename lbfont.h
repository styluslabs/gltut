#pragma once
#include "lbcommon.h"

typedef struct {
  char* ttfbuffer;
  struct stbtt_fontinfo* fontinfo;
  float EM_per_unit;
  int num_glyphs;
  LBGlyph* glyphs[128];
  GLuint shader;
} LBFont;

#include "stb_truetype.h"

void initLBGlyphPts(LBGlyph* g, float scale, int n_points, stbtt_vertex* points);
void freeLBGlyph(LBGlyph* g);
LBFont* createLBFont(const char* filename);
void freeLBFont(LBFont* font);
