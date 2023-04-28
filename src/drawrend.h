#ifndef CGL_DRAWREND_H
#define CGL_DRAWREND_H

#include "CGL/CGL.h"
#include "CGL/renderer.h"
#include "CGL/color.h"
#include <vector>
#include <cstring>
#include "GLFW/glfw3.h"
#include "svg.h"

#include "rasterizer.h"

#include <ft2build.h>
#include FT_FREETYPE_H

namespace CGL {

class DrawRend : public Renderer {
 public:
  DrawRend(std::vector<SVG*> svgs_, std::vector<FT_Face> font_faces);

  ~DrawRend( void );

  // inherited Renderer interface functions
  void init();
  void render();
  void resize( size_t w, size_t h );
  std::string name() { return "Draw"; }
  std::string info();
  void cursor_event( float x, float y );
  void scroll_event( float offset_x, float offset_y );
  void mouse_event( int key, int event, unsigned char mods );
  void keyboard_event( int key, int event, unsigned char mods );

  void set_gl(bool gl_) { gl = gl_; }

  // write current pixel buffer to disk
  void write_screenshot();

  // write only framebuffer to disk
  void write_framebuffer();

  // drawing functions
  void redraw();
  void draw_pixels();
  void draw_zoom();

  // view transform functions
  void view_init();
  void set_view(float x, float y, float span);
  void move_view(float dx, float dy, float scale);

  vector<vector<Vector2D>> drawLetter(FT_Outline *outline, float font_x, float font_y, float font_scale);
  vector<vector<Vector2D>> drawLetter(FT_Face font_face, char letter, float font_x, float font_y, float font_scale);
  FT_Outline interpolate_letter(FT_Outline *outline1, FT_Outline *outline2, float t, vector<vector<Vector2D>> pointsInContour1, vector<vector<Vector2D>> pointsInContour2);
  void drawCurve(std::vector<Vector2D> controlPoints, Color color, std::vector<Vector2D> *startingPoints, std::vector<Vector2D> *endingPoints);

  Rasterizer * software_rasterizer;

  float t = 0;
  char letter = 'A';
  char letter2 ='o';

private:
  // Global state variables for SVGs, pixels, and view transforms
  std::vector<SVG*> svgs; size_t current_svg;
  std::vector<FT_Face> font_faces;
  std::vector<Matrix3x3> svg_to_ndc;
  float view_x, view_y, view_span;

  Matrix3x3 ndc_to_screen;

  std::vector<unsigned char> framebuffer;
  size_t width, height;

  // UI state info
  float cursor_x; float cursor_y;
  bool left_clicked;
  int show_zoom;
  int sample_rate;

  PixelSampleMethod psm;
  LevelSampleMethod lsm;

  bool gl;

};

Vector2D lerp2D(Vector2D p1, Vector2D p2, float t);

} // namespace CGL

#endif // CGL_DRAWREND_H
