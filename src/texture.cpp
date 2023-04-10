#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

    float level = get_level(sp);
    int final_level = 0;

    if (sp.lsm == L_NEAREST) {
      final_level = std::lroundf(level);
      return sample_at_level(final_level, sp);
    } else if (sp.lsm == L_LINEAR) {
      int floored = std::floor(level);
      int ceiling = std::ceil(level);
      Color c1 = sample_at_level(floored, sp);
      Color c2 = sample_at_level(ceiling, sp);
      return lerp(level - (float) floored, c1, c2);
    } else {
      return sample_at_level(0, sp);
    }
  }

  Color Texture::sample_at_level(int level, const SampleParams& sp) {
    if (sp.psm == P_NEAREST) {
      return sample_nearest(sp.p_uv, level);
    } else if (sp.psm == P_LINEAR) {
      return sample_bilinear(sp.p_uv, level);
    } else {
      // return magenta for invalid level
      return Color(1, 0, 1);
    }
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    Vector2D diff_dx = (sp.p_dx_uv - sp.p_uv) * width;
    Vector2D diff_dy = (sp.p_dy_uv - sp.p_uv) * height;

    double level = std::log2(std::max(diff_dx.norm(), diff_dy.norm()));
    return (float) std::min(std::max(0.0, level), (double) mipmap.size() - 1);
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    double y = uv[1] * mip.height;
    double x = uv[0] * mip.width;
    Color texture_color = mip.get_texel((int) std::min(std::max(x, 0.0), (double) (mip.width - 1)), (int) std::min(std::max(y, 0.0), (double) (mip.height - 1)));

    return texture_color;
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    double x = uv[0] * mip.width;
    double y = uv[1] * mip.height;

    x = std::min(std::max(x, 0.0), (double) (mip.width - 1));
    y = std::min(std::max(y, 0.0), (double) (mip.height - 1));

    Vector2D u01 = Vector2D((int) ((x - 0.5) + 0.5), (int) ((y + 0.5) - 0.5)); // top left
    Vector2D u00 = Vector2D((int) ((x - 0.5) + 0.5), (int) ((y + 0.5) + 0.5)); // bottom left
    Vector2D u11 = Vector2D((int) ((x + 0.5) + 0.5), (int) ((y + 0.5) - 0.5)); // top right
    Vector2D u10 = Vector2D((int) ((x + 0.5) + 0.5), (int) ((y + 0.5) + 0.5)); // bottom right

    int u01_x = (int) std::min(std::max(u01[0], 0.0), (double) (mip.width - 1));
    int u01_y = (int) std::min(std::max(u01[1], 0.0), (double) (mip.height - 1));
    int u00_x = (int) std::min(std::max(u00[0], 0.0), (double) (mip.width - 1));
    int u00_y = (int) std::max(std::min(u00[1], (double) (mip.height - 1)), 0.0);
    int u11_x = (int) std::max(std::min(u11[0], (double) (mip.width - 1)), 0.0);
    int u11_y = (int) std::min(std::max(u11[1], 0.0), (double) (mip.height - 1));
    int u10_x = (int) std::max(std::min(u10[0], (double) (mip.width - 1)), 0.0);
    int u10_y = (int) std::max(std::min(u10[1], (double) (mip.height - 1)), 0.0);

    Color color_u00 = mip.get_texel(u00_x, u00_y);
    Color color_u10 = mip.get_texel(u10_x, u10_y);
    Color color_u01 = mip.get_texel(u01_x, u01_y);
    Color color_u11 = mip.get_texel(u11_x, u11_y);

    Vector2D e = Vector2D(x, u00[1]);
    Vector2D f_prime = Vector2D(u00[0], y);

    float alpha = (float) ((e - u00).norm() / (u10 - u00).norm());
    float beta = (float) ((f_prime - u00).norm() / (u01 - u00).norm());

    Color u1 = lerp(alpha, color_u01, color_u11);
    Color u0 = lerp(alpha, color_u00, color_u10);
    Color final = lerp(beta, u0, u1);

    return final;
  }

  /****************************************************************************/

  // Helpers
  Color lerp(float x, Color v0, Color v1) {
    return (x * v1) + ((1 - x) * v0);
    // return v0 + (x * (v1 - v0));
  }

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
