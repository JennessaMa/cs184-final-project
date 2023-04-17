#include "drawrend.h"
#include "svg.h"
#include "transforms.h"
#include "CGL/misc.h"
#include <iostream>
#include <sstream>
#include "CGL/lodepng.h"
#include "texture.h"
#include <ctime>
#include "rasterizer.h"

#define FT_ON_BIT( flag )  ( flag & 0x01 )
#define FT_ORDER_BIT( flag )  ( flag & 0x02 )

using namespace std;

namespace CGL {

struct SVG;


DrawRend::DrawRend(std::vector<SVG*> svgs_, FT_Face face)
: svgs(svgs_), current_svg(0), font_face(face)
{
}


DrawRend::~DrawRend(void) {
  svgs.clear();
  delete software_rasterizer;
}

/**
* Initialize the renderer.
* Set default parameters and initialize the viewing transforms for each tab.
*/
void DrawRend::init() {
  gl = true;

  sample_rate = 1;
  left_clicked = false;
  show_zoom = 0;

  svg_to_ndc.resize(svgs.size());
  for (int i = 0; i < svgs.size(); ++i) {
    current_svg = i;
    view_init();
  }
  current_svg = 0;
  psm = P_NEAREST;
  lsm = L_ZERO;
  
  width = height = 0;

  software_rasterizer = new RasterizerImp(psm, lsm, width, height, sample_rate);
}

/**
* Draw content.
* Simply reposts the framebuffer and the zoom window, if applicable.
*/
void DrawRend::render() {
  draw_pixels();
  if (show_zoom)
    draw_zoom();
}

/**
 * Respond to buffer resize.
 * Resizes the buffers and resets the
 * normalized device coords -> screen coords transform.
 * \param w The new width of the context
 * \param h The new height of the context
 */
void DrawRend::resize(size_t w, size_t h) {
  width = w; height = h;

  framebuffer.resize(3 * w * h);

  float scale = min(width, height);
  ndc_to_screen(0, 0) = scale; ndc_to_screen(0, 2) = (width - scale) / 2;
  ndc_to_screen(1, 1) = scale; ndc_to_screen(1, 2) = (height - scale) / 2;
  
  software_rasterizer->set_framebuffer_target(framebuffer.data(), width, height);

  redraw();
}

/**
 * Return a brief description of the renderer.
 * Displays current buffer resolution, sampling method, sampling rate.
 */
static const string level_strings[] = { "level zero", "nearest level", "bilinear level interpolation" };
static const string pixel_strings[] = { "nearest pixel", "bilinear pixel interpolation" };
std::string DrawRend::info() {
  stringstream ss;
  stringstream sample_method;
  sample_method << level_strings[lsm] << ", " << pixel_strings[psm];
  ss << "Resolution " << width << " x " << height << ". ";
  ss << "Using " << sample_method.str() << " sampling. ";
  ss << "Supersample rate " << sample_rate << " per pixel. ";
  return ss.str();
}

/**
 * Respond to cursor events.
 * The viewer itself does not really care about the cursor but it will take
 * the GLFW cursor events and forward the ones that matter to  the renderer.
 * The arguments are defined in screen space coordinates ( (0,0) at top
 * left corner of the window and (w,h) at the bottom right corner.
 * \param x the x coordinate of the cursor
 * \param y the y coordinate of the cursor
 */
void DrawRend::cursor_event(float x, float y) {
  // translate when left mouse button is held down
  if (left_clicked) {
    float dx = (x - cursor_x) / width * svgs[current_svg]->width;
    float dy = (y - cursor_y) / height * svgs[current_svg]->height;
    move_view(dx, dy, 1);
    redraw();
  }

  // register new cursor location
  cursor_x = x;
  cursor_y = y;
}

/**
 * Respond to zoom event.
 * Like cursor events, the viewer itself does not care about the mouse wheel
 * either, but it will take the GLFW wheel events and forward them directly
 * to the renderer.
 * \param offset_x Scroll offset in x direction
 * \param offset_y Scroll offset in y direction
 */
void DrawRend::scroll_event(float offset_x, float offset_y) {
  if (offset_x || offset_y) {
    float scale = 1 + 0.05 * (offset_x + offset_y);
    scale = std::min(1.5f, std::max(0.5f, scale));
    move_view(0, 0, scale);
    redraw();
  }
}

/**
 * Respond to mouse click event.
 * The viewer will always forward mouse click events to the renderer.
 * \param key The key that spawned the event. The mapping between the
 *        key values and the mouse buttons are given by the macros defined
 *        at the top of this file.
 * \param event The type of event. Possible values are 0, 1 and 2, which
 *        corresponds to the events defined in macros.
 * \param mods if any modifier keys are held down at the time of the event
 *        modifiers are defined in macros.
 */
void DrawRend::mouse_event(int key, int event, unsigned char mods) {
  if (key == MOUSE_LEFT) {
    if (event == EVENT_PRESS)
      left_clicked = true;
    if (event == EVENT_RELEASE)
      left_clicked = false;
  }
}

/**
 * Respond to keyboard event.
 * The viewer will always forward mouse key events to the renderer.
 * \param key The key that spawned the event. ASCII numbers are used for
 *        letter characters. Non-letter keys are selectively supported
 *        and are defined in macros.
 * \param event The type of event. Possible values are 0, 1 and 2, which
 *        corresponds to the events defined in macros.
 * \param mods if any modifier keys are held down at the time of the event
 *        modifiers are defined in macros.
 */
void DrawRend::keyboard_event(int key, int event, unsigned char mods) {
  if (event != EVENT_PRESS)
    return;

  // tab through the loaded files
  if (key >= '1' && key <= '9' && key - '1' < svgs.size()) {
    current_svg = key - '1';
    redraw();
    return;
  }

  switch (key) {

    // reset view transformation
  case ' ':
    view_init();
    redraw();
    break;

    // set the sampling rate to 1, 4, 9, or 16
  case '=':
    if (sample_rate < 16) {
      sample_rate = (int)(sqrt(sample_rate) + 1) * (sqrt(sample_rate) + 1);
      software_rasterizer->set_sample_rate(sample_rate);
      redraw();
    }
    break;
  case '-':
    if (sample_rate > 1) {
      sample_rate = (int)(sqrt(sample_rate) - 1) * (sqrt(sample_rate) - 1);
      software_rasterizer->set_sample_rate(sample_rate);
      redraw();
    }
    break;

    // save the current buffer to disk
  case 'S':
    write_screenshot();
    break;

    // toggle pixel sampling scheme
  case 'P':
    psm = (PixelSampleMethod)((psm + 1) % 2);
    software_rasterizer->set_psm(psm);
    redraw();
    break;
    // toggle level sampling scheme
  case 'L':
    lsm = (LevelSampleMethod)((lsm + 1) % 3);
    software_rasterizer->set_lsm(lsm);
    redraw();
    break;

    // toggle zoom
  case 'Z':
    show_zoom = (show_zoom + 1) % 2;
    break;

  default:
    return;
  }
}

/**
 * Writes the contents of the pixel buffer to disk as a .png file.
 * The image filename contains the month, date, hour, minute, and second
 * to make sure it is unique and identifiable.
 */
void DrawRend::write_screenshot() {
  redraw();
  if (show_zoom) draw_zoom();

  vector<unsigned char> windowPixels(4 * width * height);
  glReadPixels(0, 0,
    width,
    height,
    GL_RGBA,
    GL_UNSIGNED_BYTE,
    &windowPixels[0]);

  vector<unsigned char> flippedPixels(4 * width * height);
  for (int row = 0; row < height; ++row)
    memcpy(&flippedPixels[row * width * 4], &windowPixels[(height - row - 1) * width * 4], 4 * width);

  time_t t = time(nullptr);
  tm* lt = localtime(&t);
  stringstream ss;
  ss << "screenshot_" << lt->tm_mon + 1 << "-" << lt->tm_mday << "_"
    << lt->tm_hour << "-" << lt->tm_min << "-" << lt->tm_sec << ".png";
  string file = ss.str();
  cout << "Writing file " << file << "...";
  if (lodepng::encode(file, flippedPixels, width, height))
    cerr << "Could not be written" << endl;
  else
    cout << "Success!" << endl;
}

/**
 * Writes the contents of the framebuffer to disk as a .png file.
 *
 */
void DrawRend::write_framebuffer() {
  // lodepng expects alpha channel, so we will just make a new vector with
  // alpha included

  std::vector<unsigned char> export_data;

  export_data.reserve(width * height * 4);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      for (int k = 0; k < 3; ++k) {
        export_data.push_back(framebuffer[3 * (y * width + x) + k]);
      }
      export_data.push_back(255); // Opaque alpha
    }
  }

  if (lodepng::encode("test.png", export_data.data(), width, height))
    cerr << "Could not write framebuffer" << endl;
  else
    cerr << "Succesfully wrote framebuffer" << endl;
}

std::vector<Vector2D> evaluateStep(std::vector<Vector2D> const &points, float t)
{
  std::vector<Vector2D> result;
  for (int i = 0; i < points.size() - 1; i += 1) {
    Vector2D p_prime = lerp2D(points[i], points[i + 1], t);
    result.push_back(p_prime);
  }

  return result;
}

Vector2D lerp2D(Vector2D p1, Vector2D p2, float t) {
  return (1 - t) * p1 + (t * p2);
}

void DrawRend::drawCurve(std::vector<Vector2D> controlPoints)
{
  float lastX;
  float lastY;

  for (float t = 0.0; t <= 1.0f; t += 0.005f)
  {
    std::vector<Vector2D> curControlPoints = controlPoints;
    // run de casteljau until we reach 1 control point (n_points - 1 iterations)
    for (int i = 0; i < controlPoints.size() - 1; i++) {
      curControlPoints = evaluateStep(curControlPoints, t);
    }
    // curControlPoints will just have 1 point, which is the next point on the curve
    float curX = curControlPoints[0].x;
    float curY = curControlPoints[0].y;

    // draw control point
    if (t != 0.0) {
      // draw line from (lastX, lastY) to (curX, curY)
      software_rasterizer->rasterize_line(lastX, lastY, curX, curY, Color(0,0,0));
    }

    lastX = curControlPoints[0].x;
    lastY = curControlPoints[0].y;
  }
}

Vector2D offsetPoint(FT_Vector* point, Vector2D top_left, float x_offset, float y_offset) {
  float p1x = top_left.x + x_offset + (float) point->x;
  float p1y = top_left.y + y_offset + (512.0 - (float) point->y);
  return Vector2D(p1x, p1y);
}
    /**
 * Draws the current SVG tab to the screen. Also draws a
 * border around the SVG canvas. Resolves the supersample buffers
 * into the framebuffer before posting the framebuffer pixels to the screen.
 */
void DrawRend::redraw() {
  software_rasterizer->clear_buffers();

  SVG& svg = *svgs[current_svg];
  //  svg.draw(software_rasterizer, ndc_to_screen * svg_to_ndc[current_svg]);

  // draw canvas outline
  Vector2D top_left = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(0, 0)); top_left.x--; top_left.y++;
  Vector2D top_right = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(svg.width, 0)); top_right.x++; top_right.y++;
  Vector2D bottom_left = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(0, svg.height)); bottom_left.x--; bottom_left.y--;
  Vector2D bottom_right = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(svg.width, svg.height)); bottom_right.x++; bottom_right.y--;

  // Draw the font
  FT_Set_Char_Size(font_face, 0, 16 * 16, 300, 300);
  auto error = FT_Load_Glyph(font_face, FT_Get_Char_Index(font_face, 'a'), FT_LOAD_DEFAULT);
  FT_Outline *outline = &font_face->glyph->outline;

  // cout << outline->n_points << endl;
  float min_x = (float) outline->points[0].x;
  float max_x = (float) outline->points[0].x;
  float min_y = (float) outline->points[0].y;
  float max_y = (float) outline->points[0].y;
  for (int i = 1; i < outline->n_points; i += 1) {
      min_x = std::min(min_x, (float) outline->points[i].x);
      max_x = std::max(max_x, (float) outline->points[i].x);
      min_y = std::min(min_y, (float) outline->points[i].y);
      max_y = std::max(max_y, (float) outline->points[i].y);
  }

  float x_offset = ((top_right.x - top_left.x) - (max_x - min_x)) / 2.0;
  float y_offset = ((bottom_right.y - top_right.y) - (max_y - min_y)) / 2.0;

  for (int i = 0; i < outline->n_contours; i++) {
    int start_contour_index;
    if (i == 0) {
      start_contour_index = 0;
    }
    else {
      start_contour_index = outline->contours[i-1] + 1;
    }

    cout << "contour: " << start_contour_index << " to " << outline->contours[i] << endl;



  }

  std::vector<Vector2D> offsetPoints;
  std::vector<char> tags;
  for (int i = 0; i < outline->n_points; i += 1) {
    offsetPoints.push_back(offsetPoint(&outline->points[i], top_left, x_offset, y_offset));
    tags.push_back(outline->tags[i]);
  }
  offsetPoints.push_back(offsetPoint(&outline->points[0], top_left, x_offset, y_offset));
  tags.push_back(outline->tags[0]);



  for (int i = 0; i < offsetPoints.size(); i += 1) {
    // for current i, what is the start point and end point defined by the contour
    for (int i = 0; i < outline->n_contours; i++) {

    }

    // draw point at current control point
    // cout << outline->points[i].x << " " << outline->points[i].y <<  endl;
    Vector2D pt1 = offsetPoints[i];
    float p1x = pt1.x;
    float p1y = pt1.y;
    software_rasterizer->rasterize_point(p1x, p1y, Color(0, 0, 0));

    // draw line between current point and next point
    if (i < offsetPoints.size() - 1) { // 1 1
      Vector2D pt2 = offsetPoints[i + 1];
      float p2x = pt2.x;
      float p2y = pt2.y;

      if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
      FT_ON_BIT(tags[i + 1]) == FT_CURVE_TAG_ON) {
        software_rasterizer->rasterize_line(p1x, p1y, p2x, p2y, Color(0,0,0));
      }
    }

    if (i < offsetPoints.size() - 2) { // 1 0 (conic) 1
      if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
          (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 0) and
          FT_ON_BIT(tags[i + 2]) == FT_CURVE_TAG_ON) {
        // conic section
        cout << "detected conic section with 3 control points" << endl;
        std::vector<Vector2D> controlPoints;
        controlPoints.push_back(offsetPoints[i]);
        controlPoints.push_back(offsetPoints[i + 1]);
        controlPoints.push_back(offsetPoints[i + 2]);
        drawCurve(controlPoints);
      }
    }

    if (i < offsetPoints.size() - 3) { // 1 0 (cubic) 0 (cubic) 1
      if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
          (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 1) and
          (FT_ON_BIT(tags[i + 2]) == 0 and FT_ORDER_BIT(tags[i + 2]) == 1) and
          FT_ON_BIT(tags[i + 3]) == FT_CURVE_TAG_ON) {
        // conic section
        cout << "detected cubic section with 4 control points" << endl;
        std::vector<Vector2D> controlPoints;
        controlPoints.push_back(offsetPoints[i]);
        controlPoints.push_back(offsetPoints[i + 1]);
        controlPoints.push_back(offsetPoints[i + 2]);
        controlPoints.push_back(offsetPoints[i + 3]);
        drawCurve(controlPoints);
      }
    }

    if (i < offsetPoints.size() - 3) { // 1 0 (conic) 0 (conic) 1
      if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
          (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 0) and
          (FT_ON_BIT(tags[i + 2]) == 0 and FT_ORDER_BIT(tags[i + 2]) == 0) and
          FT_ON_BIT(tags[i + 3]) == FT_CURVE_TAG_ON) {
        // conic section
        cout << "detected conic off with 4 control points" << endl;
        Vector2D virtualPoint = (offsetPoints[i + 1] + offsetPoints[i + 2]) / 2;

        std::vector<Vector2D> controlPoints;
        controlPoints.push_back(offsetPoints[i]);
        controlPoints.push_back(offsetPoints[i + 1]);
        controlPoints.push_back(virtualPoint);
        drawCurve(controlPoints);

        controlPoints.clear();
        controlPoints.push_back(virtualPoint);
        controlPoints.push_back(offsetPoints[i + 2]);
        controlPoints.push_back(offsetPoints[i + 3]);
        drawCurve(controlPoints);
      }
    }

    cout << "i: " << i << " cur tag: " << FT_ON_BIT(tags[i]) << " conic or cubic: " << FT_ORDER_BIT(tags[i]) << endl;
  }

  


  // draw canvas outline
  software_rasterizer->rasterize_line(top_left.x, top_left.y, top_right.x, top_right.y, Color::Black);
  software_rasterizer->rasterize_line(top_left.x, top_left.y, bottom_left.x, bottom_left.y, Color::Black);
  software_rasterizer->rasterize_line(bottom_right.x, bottom_right.y, top_right.x, top_right.y, Color::Black);
  software_rasterizer->rasterize_line(bottom_right.x, bottom_right.y, bottom_left.x, bottom_left.y, Color::Black);

  software_rasterizer->resolve_to_framebuffer();
  if (gl)
    draw_pixels();
}

/**
 * OpenGL boilerplate to put an array of RGBA pixels on the screen.
 */
void DrawRend::draw_pixels() {
  const unsigned char* pixels = &framebuffer[0];
  // copy pixels to the screen
  glPushAttrib(GL_VIEWPORT_BIT);
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, width, 0, height, 0, 0);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(-1, 1, 0);

  glRasterPos2f(0, 0);
  glPixelZoom(1.0, -1.0);
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  glPixelZoom(1.0, 1.0);

  glPopAttrib();
  glMatrixMode(GL_PROJECTION); glPopMatrix();
  glMatrixMode(GL_MODELVIEW); glPopMatrix();
}

/**
 * Reads off the pixels that should be in the zoom window, and
 * generates a pixel array with the zoomed view.
 */
void DrawRend::draw_zoom() {

  // size (in pixels) of region of interest
  size_t regionSize = 32;

  // relative size of zoom window
  size_t zoomFactor = 16;

  // compute zoom factor---the zoom window should never cover
  // more than 40% of the framebuffer, horizontally or vertically
  size_t bufferSize = min(width, height);
  if (regionSize * zoomFactor > bufferSize * 0.4) {
    zoomFactor = (bufferSize * 0.4) / regionSize;
  }
  size_t zoomSize = regionSize * zoomFactor;

  // adjust the cursor coordinates so that the region of
  // interest never goes outside the bounds of the framebuffer
  size_t cX = max(regionSize / 2, min(width - regionSize / 2 - 1, (size_t)cursor_x));
  size_t cY = max(regionSize / 2, min(height - regionSize / 2 - 1, height - (size_t)cursor_y));

  // grab pixels from the region of interest
  vector<unsigned char> windowPixels(3 * regionSize * regionSize);
  glReadPixels(cX - regionSize / 2,
    cY - regionSize / 2 + 1, // meh
    regionSize,
    regionSize,
    GL_RGB,
    GL_UNSIGNED_BYTE,
    &windowPixels[0]);

  // upsample by the zoom factor, highlighting pixel boundaries
  vector<unsigned char> zoomPixels(3 * zoomSize * zoomSize);
  unsigned char* wp = &windowPixels[0];
  // outer loop over pixels in region of interest
  for (int y = 0; y < regionSize; y++) {
    int y0 = y * zoomFactor;
    for (int x = 0; x < regionSize; x++) {
      int x0 = x * zoomFactor;
      unsigned char* zp = &zoomPixels[(x0 + y0 * zoomSize) * 3];
      // inner loop over upsampled block
      for (int j = 0; j < zoomFactor; j++) {
        for (int i = 0; i < zoomFactor; i++) {
          for (int k = 0; k < 3; k++) {
            // highlight pixel boundaries
            if (i == 0 || j == 0) {
              const float s = .3;
              zp[k] = (int)((1. - 2. * s) * wp[k] + s * 255.);
            }
            else {
              zp[k] = wp[k];
            }
          }
          zp += 3;
        }
        zp += 3 * (zoomSize - zoomFactor);
      }
      wp += 3;
    }
  }

  // copy pixels to the screen using OpenGL
  glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity(); glOrtho(0, width, 0, height, 0.01, 1000.);
  glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity(); glTranslated(0., 0., -1.);

  glRasterPos2i(width - zoomSize, height - zoomSize);
  glDrawPixels(zoomSize, zoomSize, GL_RGB, GL_UNSIGNED_BYTE, &zoomPixels[0]);
  glMatrixMode(GL_PROJECTION); glPopMatrix();
  glMatrixMode(GL_MODELVIEW); glPopMatrix();

}

/**
 * Initializes the default viewport to center and reasonably zoom the SVG
 * with a bit of margin.
 */
void DrawRend::view_init() {
  float w = svgs[current_svg]->width, h = svgs[current_svg]->height;
  set_view(w / 2, h / 2, 1.2 * std::max(w, h) / 2);
}

/**
 * Sets the viewing transform matrix corresponding to a view centered at
 * (x,y) in SVG space, extending 'span' units in all four directions.
 * This transform maps to 'normalized device coordinates' (ndc), where the window
 * corresponds to the [0,1]^2 rectangle.
 */
void DrawRend::set_view(float x, float y, float span) {
  svg_to_ndc[current_svg] = Matrix3x3(1, 0, -x + span, 0, 1, -y + span, 0, 0, 2 * span);
}

/**
 * Recovers the previous viewing center and span from the viewing matrix,
 * then shifts and zooms the viewing window by setting a new view matrix.
 */
void DrawRend::move_view(float dx, float dy, float zoom) {
  Matrix3x3& m = svg_to_ndc[current_svg];
  float span = m(2, 2) / 2.;
  float x = span - m(0, 2), y = span - m(1, 2);
  set_view(x - dx, y - dy, span * zoom);
}

}
