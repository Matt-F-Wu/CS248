#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Extra credit: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation
  
  // return magenta for invalid level
  return Color(1,0,1,1);

}

float interpolateColorChannel(float v1, float v2, float p1, float p2, float p){
  if(p1 == p2) return v1;
  return (v2 - v1) * (p - p1) / (p2 - p1) + v1;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // return magenta for invalid level
  if( level < 0 || level >= tex.mipmap.size()) return Color(1,0,1,1);

  // Task 4: Implement bilinear filtering
  MipLevel& mip = tex.mipmap[level];
  float x = u * mip.width;
  float y = v * mip.height;

  //if (x < 0 || x >= mip.width) return Color(1.0, 1.0, 1.0, 0);
  //if (y < 0 || y >= mip.height) return Color(1.0, 1.0, 1.0, 0);

  //Get the 4 bounding colors

  float xf = x - (long)x;
  float yf = y - (long)y;

  int l = (int)floor(x);
  int r = (int)ceil(x);
  if(xf <= 0.5){
    l -= 1;
    r -= 1;
  }

  int t = (int)floor(y);
  int b = (int)ceil(y);
  if(yf <= 0.5){
    t -= 1;
    b -= 1;
  }

  float inv255 = 1.0/255;

  Color tl = Color(0.0, 0.0, 0.0, 0.0);
  Color tr = Color(0.0, 0.0, 0.0, 0.0);
  Color bl = Color(0.0, 0.0, 0.0, 0.0);
  Color br = Color(0.0, 0.0, 0.0, 0.0);

  if(t >=0 && l >=0 ) uint8_to_float(&tl.r, &mip.texels[(t*mip.width + l) * 4]);
  //cout<<tl.g<<" vs "<<mip.texels[t*mip.width + l + 1] * inv255<<endl;
  if(t >=0 && r < mip.width ) uint8_to_float(&tr.r, &mip.texels[(t*mip.width + r) * 4]);
  if(b < mip.height && l >=0 ) uint8_to_float(&bl.r, &mip.texels[(b*mip.width + l) * 4]);
  if(b < mip.height && r < mip.width ) uint8_to_float(&br.r, &mip.texels[(b*mip.width + r) * 4]);

  // Start Bilinear Algorithm
  Color top_c = Color(interpolateColorChannel(tl.r, tr.r, l+0.5, r+0.5, x), interpolateColorChannel(tl.g, tr.g, l+0.5, r+0.5, x), interpolateColorChannel(tl.b, tr.b, l+0.5, r+0.5, x), interpolateColorChannel(tl.a, tr.a, l+0.5, r+0.5, x));

  Color bottom_c = Color(interpolateColorChannel(bl.r, br.r, l+0.5, r+0.5, x), interpolateColorChannel(bl.g, br.g, l+0.5, r+0.5, x), interpolateColorChannel(bl.b, br.b, l+0.5, r+0.5, x), interpolateColorChannel(bl.a, br.a, l+0.5, r+0.5, x));

  return Color(interpolateColorChannel(top_c.r, bottom_c.r, t+0.5, b+0.5, y), interpolateColorChannel(top_c.g, bottom_c.g, t+0.5, b+0.5, y), interpolateColorChannel(top_c.b, bottom_c.b, t+0.5, b+0.5, y), interpolateColorChannel(top_c.a, bottom_c.a, t+0.5, b+0.5, y));

  // return magenta for invalid level
  // return Color(1,0,1,1);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Extra credit: Implement trilinear filtering
  float L = max(u_scale, v_scale);

  float level = log2(L);
  //The base level
  int level_b = max((int) level, 0);
  //The next level
  int level_n = max((int)level + 1, 0);

  level = max(level, 0.0f);

  //cout<<level_b<<" "<<level_n<<" actual level is: "<<level<<endl;
  Color c1 = sample_bilinear(tex, u, v, level_b);

  Color c2 = sample_bilinear(tex, u, v, level_n);

  // Interpolate c1 and c2
  return Color(interpolateColorChannel(c1.r, c2.r, level_b, level_n, level), interpolateColorChannel(c1.g, c2.g, level_b, level_n, level), interpolateColorChannel(c1.b, c2.b, level_b, level_n, level), interpolateColorChannel(c1.a, c2.a, level_b, level_n, level));

}

} // namespace CS248
