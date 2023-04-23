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
/* For geometry operations */
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>

/* For WKT read/write */
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/linearref/LengthIndexedLine.h>

using namespace geos::geom;
using namespace geos::linearref;
using namespace geos::io;

#define FT_ON_BIT( flag )  ( flag & 0x01 )
#define FT_ORDER_BIT( flag )  ( flag & 0x02 )

using namespace std;

namespace CGL {

struct SVG;


DrawRend::DrawRend(std::vector<SVG*> svgs_, std::vector<FT_Face> font_faces_)
: svgs(svgs_), current_svg(0), font_faces(font_faces_)
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

bool pointInPolygon(int polyCorners, vector<Vector2D> startingPoints, vector<Vector2D> endingPoints, float x, float y) {
  bool inside = false;

  for (int i=0; i < polyCorners; i++) {
    Vector2D start = startingPoints[i];
    Vector2D end = endingPoints[i];
    if ((start.y < y && end.y >= y) || (end.y < y && start.y >= y)) {
      if (start.x + (y - start.y) / (end.y - start.y) * (end.x - start.x) > x) {
        inside = !inside;
      }
    }
  }

  return inside;
}

void DrawRend::drawCurve(std::vector<Vector2D> controlPoints, Color color, std::vector<Vector2D> *startingPoints, std::vector<Vector2D> *endingPoints)
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
      software_rasterizer->rasterize_line(lastX, lastY, curX, curY, color);
      startingPoints->push_back(Vector2D(lastX, lastY));
      endingPoints->push_back(Vector2D(curX, curY));
    }

    lastX = curControlPoints[0].x;
    lastY = curControlPoints[0].y;
  }
}

FT_Outline DrawRend::interpolate_letter(FT_Outline *outline1, FT_Outline *outline2, float t, vector<vector<Vector2D>> pointsInContour1, vector<vector<Vector2D>> pointsInContour2) {
  FT_Outline res;
  int m = 100;
  res.n_contours = outline1->n_contours;
  res.n_points = 0; // FIXME: update to # control points

  vector<vector<Vector2D>> sampledPointsPerContour1; // list of m points sampled per contour
  vector<vector<Vector2D>> sampledPointsPerContour2; // list of m points sampled per contour

  GeometryFactory::Ptr factory = GeometryFactory::create();

  for (int i = 0; i < outline1->n_contours; i += 1) {
    // sample m points on the first font
    vector<Vector2D> contour1Points = pointsInContour1[i];
    vector<Vector2D> mSampledPointsForCurContour1;

    CoordinateArraySequence *cl1 = new CoordinateArraySequence();
    // add all the contour 1 points to cl
    for (auto pt : contour1Points) {
      cl1->add(Coordinate(pt.x, pt.y));
    }

    LineString *ls1 = factory->createLineString(cl1);
    // compute contour length and spacing
    float spacing1 = ls1->getLength() / m;
    LengthIndexedLine *lin1 = new LengthIndexedLine(ls1);

    for (int j = 0; j < m; j++) {
      // get point j along the contour
      Coordinate newPoint = lin1->extractPoint(0, j*spacing1);
      mSampledPointsForCurContour1.push_back(Vector2D(newPoint.x, newPoint.y));
    }

    // sample m points on the second font
    vector<Vector2D> contour2Points = pointsInContour2[i];
    vector<Vector2D> mSampledPointsForCurContour2;
    CoordinateArraySequence *cl2 = new CoordinateArraySequence();

    // add all the contour 2 points to cl
    for (auto pt : contour2Points) {
      cl2->add(Coordinate(pt.x, pt.y));
    }

    LineString *ls2 = factory->createLineString(cl2);
    // compute contour length and spacing
    float spacing2 = ls2->getLength() / m;
    LengthIndexedLine *lin2 = new LengthIndexedLine(ls2);

    for (int j = 0; j < m; j++) {
      // get point j along the contour 2
      Coordinate newPoint = lin2->extractPoint(0, j*spacing2);
      mSampledPointsForCurContour2.push_back(Vector2D(newPoint.x, newPoint.y));
    }

    // TODO: find index in mSampledPointsForCurContour1 and mSampledPointsForCurContour2 that would be the starting point
    int contour1Offset = 0;
    int contour2Offset = 0;

    // lerp the 2 sets of m sampled points
    vector<Vector2D> lerpedPoints;
    for (int j = 0; j < m; j++) {
      // TODO: offset j by indices % 100 (starting point)
      Vector2D interpolatedPoint = lerp2D(mSampledPointsForCurContour1[j], mSampledPointsForCurContour2[j], t);
      lerpedPoints.push_back(interpolatedPoint);
    }

    // TODO: fit bezier curves to 100 lerped points (find control points)
  }

  return *outline1;
}

vector<vector<Vector2D>> DrawRend::drawLetter(FT_Outline *outline, float font_x, float font_y, float font_scale) {
  SVG& svg = *svgs[current_svg];
  //  svg.draw(software_rasterizer, ndc_to_screen * svg_to_ndc[current_svg]);

  // draw canvas outline
  Vector2D top_left = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(0, 0)); top_left.x--; top_left.y++;
  Vector2D top_right = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(svg.width, 0)); top_right.x++; top_right.y++;
  Vector2D bottom_left = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(0, svg.height)); bottom_left.x--; bottom_left.y--;
  Vector2D bottom_right = ndc_to_screen * svg_to_ndc[current_svg] * (Vector2D(svg.width, svg.height)); bottom_right.x++; bottom_right.y--;


  vector<Color> contourColors;
  contourColors.push_back(Color(1, 0, 0)); //red
  contourColors.push_back(Color(0, 1, 0)); //green
  contourColors.push_back(Color(0, 0, 1)); //blue
  contourColors.push_back(Color(1, 1, 0)); //yellow
  contourColors.push_back(Color(1, 1, 1)); //black
  contourColors.push_back(Color(1, 0, 1)); //pink

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

  float canvas_width = top_right.x - top_left.x;
  float canvas_height = bottom_right.y - top_right.y;

  vector<Vector2D> shiftedPoints;
  for (int i = 0; i < outline->n_points; i += 1) {
    int x_width = max_x - min_x;
    int y_height = max_y - min_y;
    int shared_size = max(x_width, y_height);
    double shiftX = (outline->points[i].x - min_x) / shared_size;
    double shiftY = (outline->points[i].y - min_y) / shared_size;
    shiftedPoints.push_back(Vector2D(font_x * canvas_width + top_left.x + shiftX * canvas_width * font_scale,
                                     font_y * canvas_height + top_left.y + (canvas_height - shiftY * canvas_height) * font_scale));
  }
  vector<Vector2D> startingPoints;
  vector<Vector2D> endingPoints;

  vector<int> numPointsPerContour;

  for (int j = 0; j < outline->n_contours; j++) {
    int numStartPointsDrawnSoFar = startingPoints.size();
    int start_contour_index;
    if (j == 0) {
      start_contour_index = 0;
    }
    else {
      start_contour_index = outline->contours[j-1] + 1;
    }
    int end_contour_index = outline->contours[j];

    // cout << "contour #" << j << ": " << start_contour_index << " to " << end_contour_index << endl;
    int contour_length = end_contour_index - start_contour_index;

    std::vector<Vector2D> curPoints; // length of this is contour_length
    std::vector<char> tags; // length of this is contour_length
    for (int i = start_contour_index; i <= end_contour_index; i += 1) {
      curPoints.push_back(shiftedPoints[i]);
      tags.push_back(outline->tags[i]);
    }
    curPoints.push_back(shiftedPoints[start_contour_index]);
    tags.push_back(outline->tags[start_contour_index]);

    for (int i = 0; i < curPoints.size(); i += 1) {
      // draw point at current control point
      // cout << outline->points[i].x << " " << outline->points[i].y <<  endl;
      Vector2D pt1 = curPoints[i];
      float p1x = pt1.x;
      float p1y = pt1.y;
      software_rasterizer->rasterize_point(p1x, p1y, Color(0, 0, 0));

      // draw line between current point and next point
      if (i < curPoints.size() - 1) { // 1 1
        Vector2D pt2 = curPoints[i + 1];
        float p2x = pt2.x;
        float p2y = pt2.y;

        if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
            FT_ON_BIT(tags[i + 1]) == FT_CURVE_TAG_ON) {
          software_rasterizer->rasterize_line(p1x, p1y, p2x, p2y, contourColors[j]);
          startingPoints.push_back(Vector2D(p1x, p1y));
          endingPoints.push_back(Vector2D(p2x, p2y));
        }
      }

      if (i < curPoints.size() - 2) { // 1 0 (conic) 1
        if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
            (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 0) and
            FT_ON_BIT(tags[i + 2]) == FT_CURVE_TAG_ON) {
          // conic section
          // cout << "detected conic section with 3 control points" << endl;
          std::vector<Vector2D> controlPoints;
          controlPoints.push_back(curPoints[i]);
          controlPoints.push_back(curPoints[i + 1]);
          controlPoints.push_back(curPoints[i + 2]);
          drawCurve(controlPoints, contourColors[j], &startingPoints, &endingPoints);
        }
      }

      if (i < curPoints.size() - 3) { // 1 0 (cubic) 0 (cubic) 1
        if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
            (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 1) and
            (FT_ON_BIT(tags[i + 2]) == 0 and FT_ORDER_BIT(tags[i + 2]) == 1) and
            FT_ON_BIT(tags[i + 3]) == FT_CURVE_TAG_ON) {
          // conic section
          // cout << "detected cubic section with 4 control points" << endl;
          std::vector<Vector2D> controlPoints;
          controlPoints.push_back(curPoints[i]);
          controlPoints.push_back(curPoints[i + 1]);
          controlPoints.push_back(curPoints[i + 2]);
          controlPoints.push_back(curPoints[i + 3]);
          drawCurve(controlPoints, contourColors[j], &startingPoints, &endingPoints);
        }
      }

      if (i < curPoints.size() - 4) { // 1 0 (conic) 0 (conic) 0 (conic) 1
        if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
            (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 0) and
            (FT_ON_BIT(tags[i + 2]) == 0 and FT_ORDER_BIT(tags[i + 2]) == 0) and
            (FT_ON_BIT(tags[i + 3]) == 0 and FT_ORDER_BIT(tags[i + 3]) == 0) and
            FT_ON_BIT(tags[i + 4]) == FT_CURVE_TAG_ON) {
          // conic section
          // cout << "detected cubic section with 4 control points" << endl;
          std::vector<Vector2D> controlPoints;
          controlPoints.push_back(curPoints[i]);
          controlPoints.push_back(curPoints[i + 1]);
          controlPoints.push_back(curPoints[i + 2]);
          controlPoints.push_back(curPoints[i + 3]);
          controlPoints.push_back(curPoints[i + 4]);
          drawCurve(controlPoints, contourColors[j], &startingPoints, &endingPoints);
        }
      }

      if (i < curPoints.size() - 3) { // 1 0 (conic) 0 (conic) 1
        if (FT_ON_BIT(tags[i]) == FT_CURVE_TAG_ON and
            (FT_ON_BIT(tags[i + 1]) == 0 and FT_ORDER_BIT(tags[i + 1]) == 0) and
            (FT_ON_BIT(tags[i + 2]) == 0 and FT_ORDER_BIT(tags[i + 2]) == 0) and
            FT_ON_BIT(tags[i + 3]) == FT_CURVE_TAG_ON) {
          // conic section
          // cout << "detected conic off with 4 control points" << endl;
          Vector2D virtualPoint = (curPoints[i + 1] + curPoints[i + 2]) / 2;

          std::vector<Vector2D> controlPoints;
          controlPoints.push_back(curPoints[i]);
          controlPoints.push_back(curPoints[i + 1]);
          controlPoints.push_back(virtualPoint);
          drawCurve(controlPoints, contourColors[j], &startingPoints, &endingPoints);

          controlPoints.clear();
          controlPoints.push_back(virtualPoint);
          controlPoints.push_back(curPoints[i + 2]);
          controlPoints.push_back(curPoints[i + 3]);
          drawCurve(controlPoints, contourColors[j], &startingPoints, &endingPoints);
        }
      }

      // cout << "i: " << i + start_contour_index << " cur tag: " << FT_ON_BIT(tags[i]) << " conic or cubic: " << FT_ORDER_BIT(tags[i]) << endl;
    }

    int numPointsForThisContour = startingPoints.size() - numStartPointsDrawnSoFar;
    numPointsPerContour.push_back(numPointsForThisContour);
  }

  // fill in glyph
  for (int x = (int) top_left.x; x < (int) top_right.x; x++) {
    for (int y = (int) top_left.y; y < (int) bottom_left.y; y++) {
      if (pointInPolygon(startingPoints.size(), startingPoints, endingPoints, (float) x, (float) y)) {
        software_rasterizer->rasterize_point((float) x, (float) y, Color::Black);
      }
    }
  }

  vector<vector<Vector2D>> startingPointsPerContour;

  int contour_start_i = 0;
  int contour_end_i = 0;
  for (int i = 0; i < outline->n_contours; i++) {
    contour_end_i += numPointsPerContour[i];

    // take points from contour_start_i to contour_end_i (these are the ones that are on the current contour)
    vector<Vector2D> pointsOnCurrentContour;
    for (int j = contour_start_i; j < contour_end_i; j++) {
      pointsOnCurrentContour.push_back(startingPoints[i]);
    }
    startingPointsPerContour.push_back(pointsOnCurrentContour);

    contour_start_i += numPointsPerContour[i];
  }

  return startingPointsPerContour;
}

vector<vector<Vector2D>> DrawRend::drawLetter(FT_Face font_face, char letter, float font_x, float font_y, float font_scale) {
  // Draw the font
  FT_Set_Char_Size(font_face, 0, 16 * 16, 300, 300);
  auto error = FT_Load_Glyph(font_face, FT_Get_Char_Index(font_face, letter), FT_LOAD_DEFAULT);
  FT_Outline *outline = &font_face->glyph->outline;

  return drawLetter(outline, font_x, font_y, font_scale);
}

    /**
 * Draws the current SVG tab to the screen. Also draws a
 * border around the SVG canvas. Resolves the supersample buffers
 * into the framebuffer before posting the framebuffer pixels to the screen.
 */
void DrawRend::redraw() {
  software_rasterizer->clear_buffers();

  if (font_faces.size() == 1) {
    // draw just a single font
    drawLetter(font_faces[0], 'A', 0, 0, 0.5);
  } else if (font_faces.size() == 2){
    // draw two fonts and interpolate
    char letter = 'A';
    FT_Set_Char_Size(font_faces[0], 0, 16 * 16, 300, 300);
    FT_Load_Glyph(font_faces[0], FT_Get_Char_Index(font_faces[0], letter), FT_LOAD_DEFAULT);
    FT_Outline *outline1 = &font_faces[0]->glyph->outline;
    cout << outline1->n_points << endl;

    FT_Set_Char_Size(font_faces[1], 0, 16 * 16, 300, 300);
    FT_Load_Glyph(font_faces[1], FT_Get_Char_Index(font_faces[1], letter), FT_LOAD_DEFAULT);
    FT_Outline *outline2 = &font_faces[1]->glyph->outline;
    cout << outline2->n_points << endl;

    vector<vector<Vector2D>> first_letter_points = drawLetter(font_faces[0], letter, 0, 0, 0.33);
    vector<vector<Vector2D>> second_letter_points = drawLetter(font_faces[1], letter, 0.66, 0, 0.33);
//    FT_Outline interpolated_outline = interpolate_letter(outline1, outline2, 0.5, first_letter_points, second_letter_points);
//    drawLetter(&interpolated_outline, 0.33, 0, 0.33);

  }

  // draw canvas outline
//  software_rasterizer->rasterize_line(top_left.x, top_left.y, top_right.x, top_right.y, Color::Black);
//  software_rasterizer->rasterize_line(top_left.x, top_left.y, bottom_left.x, bottom_left.y, Color::Black);
//  software_rasterizer->rasterize_line(bottom_right.x, bottom_right.y, top_right.x, top_right.y, Color::Black);
//  software_rasterizer->rasterize_line(bottom_right.x, bottom_right.y, bottom_left.x, bottom_left.y, Color::Black);

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
