#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    int sqrt_sample_rate = (int) sqrt((double) this->sample_rate);

    for (int i = x * sqrt_sample_rate; i < (x+1) * sqrt_sample_rate; i++) {
      for (int j = y * sqrt_sample_rate; j < (y+1) * sqrt_sample_rate; j++) {
        sample_buffer[(j * width * sqrt_sample_rate) + i] = c;
      }
    }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
//    float startX = min(min(x0, x1), x2);
//    float startY = min(min(y0, y1), y2);
//    float maxX = max(max(x0, x1), x2);
//    float maxY = max(max(y0, y1), y2);
//    for (int i = startX; i < maxX; i += 1) {
//        for (int j = startY; j < maxY; j += 1) {
//
//          bool test1 = inside(x0, y0, x1, y1, i + 0.5, j + 0.5);
//          bool test2 = inside(x1, y1, x2, y2, i + 0.5, j + 0.5);
//          bool test3 = inside(x2, y2, x0, y0, i + 0.5, j + 0.5);
//          if ((test1 && test2 && test3) || (!test1 && !test2 && !test3)) {
//            fill_pixel((size_t) i, (size_t) j, color);
//          }
//        }
//    }
    // TODO: Task 2: Update to implement super-sampled rasterization
    rasterize_interpolated_color_triangle(x0,y0,color,x1,y1,color,x2,y2,color);
  }

  bool RasterizerImp::inside(float x0, float y0, float x1, float y1, float x, float y) {
    float dX = x1-x0;
    float dY = y1-y0;

    float result = -(x-x0)*dY + (y-y0)*dX;

    return result >= 0;
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    int sqrt_sample_rate = (int) sqrt((double) this->sample_rate);

    x0 = x0 * sqrt_sample_rate;
    y0 = y0 * sqrt_sample_rate;
    x1 = x1 * sqrt_sample_rate;
    y1 = y1 * sqrt_sample_rate;
    x2 = x2 * sqrt_sample_rate;
    y2 = y2 * sqrt_sample_rate;
    float startX = min(min(x0, x1), x2);
    float startY = min(min(y0, y1), y2);
    float maxX = max(max(x0, x1), x2);
    float maxY = max(max(y0, y1), y2);
    for (int i = startX; i < maxX; i += 1) {
      for (int j = startY; j < maxY; j += 1) {
        bool test1 = inside(x0, y0, x1, y1, i + 0.5, j + 0.5);
        bool test2 = inside(x1, y1, x2, y2, i + 0.5, j + 0.5);
        bool test3 = inside(x2, y2, x0, y0, i + 0.5, j + 0.5);
        if ((test1 && test2 && test3) || (!test1 && !test2 && !test3)) {
          Matrix3x3 M = Matrix3x3(x0,x1,x2,
                                  y0,y1,y2,
                                  1,1,1);
          Vector3D v = Vector3D(i,j,1);
          Vector3D bary = M.inv()*v;
          Color cur_color = bary[0] * c0 + bary[1] * c1 + bary[2] * c2;

          sample_buffer[(j * width * sqrt_sample_rate) + i] = cur_color;
        }
      }
    }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    int sqrt_sample_rate = (int) sqrt((double) this->sample_rate);

    x0 = x0 * sqrt_sample_rate;
    y0 = y0 * sqrt_sample_rate;
    x1 = x1 * sqrt_sample_rate;
    y1 = y1 * sqrt_sample_rate;
    x2 = x2 * sqrt_sample_rate;
    y2 = y2 * sqrt_sample_rate;
    float startX = min(min(x0, x1), x2);
    float startY = min(min(y0, y1), y2);
    float maxX = max(max(x0, x1), x2);
    float maxY = max(max(y0, y1), y2);
    for (int i = startX; i < maxX; i += 1) {
      for (int j = startY; j < maxY; j += 1) {
        bool test1 = inside(x0, y0, x1, y1, i + 0.5, j + 0.5);
        bool test2 = inside(x1, y1, x2, y2, i + 0.5, j + 0.5);
        bool test3 = inside(x2, y2, x0, y0, i + 0.5, j + 0.5);
        if ((test1 && test2 && test3) || (!test1 && !test2 && !test3)) {

          Matrix3x3 M = Matrix3x3(x0,x1,x2,
                                  y0,y1,y2,
                                  1,1,1);
          Vector3D v = Vector3D(i,j,1);
          Vector3D bary = M.inv()*v;
          Vector2D p_uv = bary[0] * Vector2D(u0, v0) + bary[1] * Vector2D(u1, v1) + bary[2] * Vector2D(u2, v2);

          Vector3D v_dx =  Vector3D(i + 1,j,1);
          Vector3D bary_dx = M.inv() * v_dx;
          Vector2D p_dx = bary_dx[0] * Vector2D(u0, v0) + bary_dx[1] * Vector2D(u1, v1) + bary_dx[2] * Vector2D(u2, v2);

          Vector3D v_dy =  Vector3D(i,j + 1,1);
          Vector3D bary_dy = M.inv() * v_dy;
          Vector2D p_dy = bary_dy[0] * Vector2D(u0, v0) + bary_dy[1] * Vector2D(u1, v1) + bary_dy[2] * Vector2D(u2, v2);

          SampleParams sp;
          sp.p_uv = p_uv;
          sp.p_dx_uv = p_dx;
          sp.p_dy_uv = p_dy;
          sp.lsm = lsm;
          sp.psm = psm;

          // cout << "THIS IS U0, V0: " << u0 <<", " << v0 << endl;
          sample_buffer[(j * width * sqrt_sample_rate) + i] = tex.sample(sp);
        }
      }
    }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    this->sample_buffer.resize(width * height * rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support

    int sqrt_sample_rate = (int) sqrt((double) this->sample_rate);
    int sample_buffer_width = width * sqrt_sample_rate;

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color avg = Color();
        for (int i = x * sqrt_sample_rate; i < (x+1) * sqrt_sample_rate; i++) {
          for (int j = y * sqrt_sample_rate; j < (y+1) * sqrt_sample_rate; j++) {
            avg += sample_buffer[sample_buffer_width * j + i];
          }
        }
        avg = avg * (1.0 / sample_rate);
        // TODO: do averaging of supersamples here, final result goes to rgb_framebuffer_target

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&avg.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
