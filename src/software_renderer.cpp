#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

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
  dst_uint8[0] = (uint8_t) floor( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) floor( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) floor( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) floor( 255.f * max( 0.0f, min( 1.0f, src[3])));
}
// My way of determining if normal vector is an inner normal
bool isInner(float normal_x, float normal_y, float other_edge_x, float other_edge_y){
  if(normal_x * other_edge_x + normal_y * other_edge_y >= 0){
    return true;
  }else{
    return false;
  }
}

// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int x, int y, const Color &color) {
  // Hao's optimization, I am maximizing memory efficiency here

  // Fill sample shouldn't really fill pixels directly, but filling the supersample_target
  //fill_pixel((int)floor(sx), (int)floor(sy), color);
  if (x < 0 || x >= target_w * sample_rate) return;
  if (y < 0 || y >= target_h * sample_rate) return;

  // Pixel is out of bound, return
  //if (x < 0 || x >= target_w * target_h * sample_rate * sample_rate * 4) {
  //  return;
  //}

  Color pixel_color;
  uint8_to_float(&pixel_color.r, &supersample_target[(x + y * target_w * sample_rate)*4]);

  pixel_color = ref->alpha_blending_helper(pixel_color, color);

  float_to_uint8(&supersample_target[(x + y * target_w * sample_rate)*4], &pixel_color.r);

  //cout << "super sample color: " << pixel_color.r * 255 << endl;
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= target_w) return;
	if (y < 0 || y >= target_h) return;

  if(resolving){

  	Color pixel_color;
    uint8_to_float(&pixel_color.r, &render_target[4 * (x + y * target_w)]);

  	pixel_color = ref->alpha_blending_helper(pixel_color, color);

  	float_to_uint8(&render_target[4 * (x + y * target_w)], &pixel_color.r);

    return;
  }
  /* fill sample buffer with the pixel you draw, so points and lines that are directly drawn to 
     render_target will not be lost during the resolve step (resampling supersample_target to render_target)
  */

  for(int i = 0; i < sample_rate; i++){
    for(int j = 0; j < sample_rate; j++){
      //fill corresponding sample locations
      //fill_sample(4 * ((x + y * target_w) * (sample_rate * sample_rate) + i * sample_rate + j), 0, color);
      fill_sample(x * sample_rate + i, y * sample_rate + j, color);
    }
  }

  //delete &color;
}

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = canvas_to_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  // How to only resolve triangles, but leaving points and lines alone
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  unsigned char* tmp = supersample_target;
  supersample_target = new unsigned char[target_w * target_h * sample_rate * sample_rate * 4];
  delete tmp;
  clear_target();
  memset(supersample_target, 255, 4 * target_w * target_h * sample_rate * sample_rate);
}

//set_render_target() is called whenever the user resizes the application window.
void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  clear_target();
  unsigned char* tmp = supersample_target;
  supersample_target = new unsigned char[target_w * target_h * sample_rate * sample_rate * 4];
  if(tmp != NULL){
    delete tmp;
  }
  memset(supersample_target, 255, 4 * target_w * target_h * sample_rate * sample_rate);
  //cout << "Init supersample_target: " << (int)(supersample_target[0]) << endl;
}


void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack

  // Achieves a stack push
  Matrix3x3 old_transformation = transformation;
  transformation = transformation * element->transform;

	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}

  // Esssentially a stack pop
  transformation = old_transformation;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel (Actually filling in the pixel the sample takes place in)
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= target_w) return;
  if (sy < 0 || sy >= target_h) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  fill_pixel(sx, sy, color);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Extra credit (delete the line below and implement your own)
  ref->rasterize_line_helper(x0, y0, x1, y1, target_w, target_h, color, this);

}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 1: 
  // Implement triangle rasterization (you may want to call fill_sample here)
  
  //First find the bounding box
  float max_x = ceil(max(max(x0, x1), x2));
  float max_y = ceil(max(max(y0, y1), y2));

  float min_x = floor(min(min(x0, x1), x2));
  float min_y = floor(min(min(y0, y1), y2));

  // Then find the 3 edges' inner norm vector
  float e0_nor_x = y1 - y0;
  float e0_nor_y = -(x1 - x0);

  float e1_nor_x = y2 - y1;
  float e1_nor_y = -(x2 - x1);

  float e2_nor_x = y0 - y2;
  float e2_nor_y = -(x0 - x2);
  
  if(!isInner(e0_nor_x, e0_nor_y, x2 - x0, y2 - y0)){
    e0_nor_x *= -1;
    e0_nor_y *= -1;
  }

  if(!isInner(e1_nor_x, e1_nor_y, x0 - x1, y0 - y1)){
    e1_nor_x *= -1;
    e1_nor_y *= -1;
  }

  if(!isInner(e2_nor_x, e2_nor_y, x1 - x2, y1 - y2)){
    e2_nor_x *= -1;
    e2_nor_y *= -1;
  }

  //Now we have all the inner normals, we can just verify each point in the bounding box
  float x_s, y_s;
  for(int y = min_y; y < max_y; y++){
    for(int x = min_x; x < max_x; x++){
      //The step size of super samples
      float super_sample_step_size = 1.0/(2*sample_rate);
      //Location of the top left super sample of a pixel
      x_s = x + super_sample_step_size;
      y_s = y + super_sample_step_size;

      for(int j = 0; j < sample_rate; j++){
        for(int i = 0; i < sample_rate; i++){
          if(isInner(e0_nor_x, e0_nor_y, x_s + i*super_sample_step_size*2 - x0, y_s + j*super_sample_step_size*2 - y0) &&
            isInner(e1_nor_x, e1_nor_y, x_s + i*super_sample_step_size*2 - x1, y_s + j*super_sample_step_size*2 - y1) &&
            isInner(e2_nor_x, e2_nor_y, x_s + i*super_sample_step_size*2 - x2, y_s + j*super_sample_step_size*2 - y2)){
            //The point is inside the triangle
            //cout << "Inside..." << endl;
            //fill_sample(4 * ((x + y * target_w) * (sample_rate * sample_rate) + i * sample_rate + j), 0, color);
            fill_sample(x * sample_rate + i, y * sample_rate + j, color);
          }

        }  
      }
        
    }
  }
  
  return;
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization (you may want to call fill_sample here)
  double data[9] = {1.0/(x1-x0), 0, 1.0/(x1-x0)* -1.0 * x0, 0, 1.0/(y1-y0), 1.0/(y1-y0) * -1.0 * y0, 0, 0, 1};
  //Store old transformation, so we can recover later
  Matrix3x3 old_transformation = transformation;
  transformation = Matrix3x3(data);

  //Create sampler
  Sampler2D* s = new Sampler2DImp(BILINEAR);

  float super_sample_step_size = 1.0/(2*sample_rate);
  float x_s, y_s;

  for(int y = y0; y <= y1; y++){
    for(int x = x0; x <= x1; x++){
      //Get the coordinate in texture space
      x_s = x + super_sample_step_size;
      y_s = y + super_sample_step_size;

      for(int j = 0; j < sample_rate; j++){
        for(int i = 0; i < sample_rate; i++){

          Vector2D p = transform(Vector2D( x_s + i*super_sample_step_size*2, y_s + j*super_sample_step_size*2));
          
          //Bilinear sample tex
          Color c = s->sample_bilinear(tex, p.x, p.y, 0);

          /*Trilinear sampling code*/
          // Get the top sample's mapped location
          Vector2D p_t = transform(Vector2D( x_s + i*super_sample_step_size*2, y_s + (j-1)*super_sample_step_size*2 ));
          // Get the right sample's mapped location
          Vector2D p_r = transform(Vector2D( x_s + (i+1)*super_sample_step_size*2, y_s + j*super_sample_step_size*2 ));

          float du_dx = (p_r.x - p.x)*tex.mipmap[0].width / (2*super_sample_step_size);
          //cout<<du_dx<<endl;
          float dv_dx = (p_r.y - p.y)*tex.mipmap[0].height / (2*super_sample_step_size);
          //cout<<dv_dx<<endl;
          float du_dy = (p.x - p_t.x)*tex.mipmap[0].width / (2*super_sample_step_size);
          //cout<<du_dy<<endl;
          float dv_dy = (p.y - p_t.y)*tex.mipmap[0].height / (2*super_sample_step_size); 
          //cout<<dv_dy<<"======="<<endl;
          c = s->sample_trilinear(tex, p.x, p.y, sqrt(du_dx*du_dx + dv_dx*dv_dx), sqrt(du_dy*du_dy + dv_dy*dv_dy));

          //fill_sample(4 * ((x + y * target_w) * (sample_rate * sample_rate) + i * sample_rate + j), 0, c);
          fill_sample(x * sample_rate + i, y * sample_rate + j, c);

        }
      }
    }
  }

  transformation = old_transformation;
  // Why is this causing error?
  //delete s;
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
  resolving = true;
  //Box filter on supersample_target, stride is this -> sample_rate
  for (int x = 0; x < target_w; x++)
  {
    for(int y = 0; y < target_h; y++){
      //look at all the actual sample locations, and box filter on the super samples contained in each sample

      // Don't forget to initialize these values, or will have unexpected behavior
      float r_sum = 0.0, g_sum = 0.0, b_sum = 0.0, a_sum = 0.0;
      for(int j = 0; j < sample_rate; j++){
        for(int i = 0; i < sample_rate; i++){
          r_sum += (int)supersample_target[4 * (x * sample_rate + i + (y * sample_rate + j) * target_w * sample_rate)];
          g_sum += (int)supersample_target[4 * (x * sample_rate + i + (y * sample_rate + j) * target_w * sample_rate) + 1];
          b_sum += (int)supersample_target[4 * (x * sample_rate + i + (y * sample_rate + j) * target_w * sample_rate) + 2];
          a_sum += (int)supersample_target[4 * (x * sample_rate + i + (y * sample_rate + j) * target_w * sample_rate) + 3];
        }
      }
      //get the r,g,b,a average
      r_sum /= (sample_rate*sample_rate);
      g_sum /= (sample_rate*sample_rate);
      b_sum /= (sample_rate*sample_rate);
      a_sum /= (sample_rate*sample_rate);

      float inv255 = 1.0 / 255.0;

      Color c = Color(r_sum * inv255, g_sum * inv255, b_sum * inv255, a_sum * inv255);
      
      //cout << "fill_pixel: " << c->r << "  supersample_target: " << (int)supersample_target[4 * (x * sample_rate + 0 + (y * sample_rate + 0) * target_w * sample_rate)] << endl;

      fill_pixel(x, y, c);
    }
  }

  // reset flag
  resolving = false;

  memset(supersample_target, 255, 4 * target_w * target_h * sample_rate * sample_rate);

  return;

}


} // namespace CS248
