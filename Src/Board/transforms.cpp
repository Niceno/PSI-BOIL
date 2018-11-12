#include "rect.h"
#include "shapes.h"
#include "transforms.h"
#include <cmath>

//
// Transform
// 

real
Transform::rounded( real x ) const
{
  return round( 10000*x ) / 10000;
} 

real
Transform::mapX( real x ) const
{
  return rounded( x * _scale + _deltaX );
}

real
Transform::scale( real x ) const
{
  return rounded( x * _scale );
}

void
Transform::apply( real & x, real & y ) const
{
  x = mapX( x );
  y = mapY( y );
}

void
Transform::furtherScale( real x )
{
  _scale *= x;
}

//
// TransformEPS
// 

real
TransformEPS::mapY( real y ) const
{
  return rounded( y * _scale + _deltaY );
}

void
TransformEPS::setBoundingBox( const Rect & rect )
{
  real ppmm = (720/254.0f);  
  if ( rect.height/rect.width > 27.7/19.0 ) {
    _scale = 277 * ppmm / rect.height;
  } else {
    _scale = 190 * ppmm / rect.width;
  }
  _deltaX = 0.5 * 210 * ppmm - _scale * ( rect.left + 0.5 * rect.width );
  _deltaY = 0.5 * 297 * ppmm - _scale * ( rect.top - 0.5 * rect.height );  
}

//
// TransformFIG
// 

real
TransformFIG::rounded( real x ) const
{
  return round( x );
}

real
TransformFIG::mapY( real y ) const
{
  // real ppmm = 1200/25.4;
  real ppmm = 1143/25.4f;
  return rounded( ( 297*ppmm ) - ( y * _scale + _deltaY ) );
}

int
TransformFIG::mapWidth( real width ) const
{
  // FIG width unit is 1/160 inch
  // Postscript points are 1/72 inch
  int result =  static_cast<int>( round( 160 * (width / 72.0 ) ) );
  return result?result:1;
}

void
TransformFIG::setBoundingBox( const Rect & rect )
{
  // real ppmm = (1200/25.4);
  real ppmm = (1143/25.4f);
  if ( rect.height/rect.width > 27.7/19.0 ) {
    _scale = 277 * ppmm / rect.height;
  } else {
    _scale = ( 190 * ppmm ) / rect.width;
  }
  _deltaX = 0.5 * 210 * ppmm - _scale * ( rect.left + 0.5 * rect.width );
  _deltaY = 0.5 * 297 * ppmm - _scale * ( rect.top - 0.5 * rect.height );  
}

void
TransformFIG::setDepthRange( const std::vector<Form*> & shapes )
{
  maxDepth = 0;
  minDepth = std::numeric_limits<unsigned int>::max();
  std::vector<Form*>::const_iterator i = shapes.begin();
  std::vector<Form*>::const_iterator end = shapes.end();
  while ( i != end ) {
    if ( (*i)->depth > maxDepth ) maxDepth = (*i)->depth;
    if ( (*i)->depth < minDepth ) minDepth = (*i)->depth;
    ++i;
  }
}

unsigned int
TransformFIG::mapDepth( unsigned int depth ) const
{
  if ( maxDepth - minDepth > 998 ) {
    real range = maxDepth - minDepth;
    return static_cast<unsigned int>( 1 + round( ( ( depth - minDepth ) / range ) * 998 ) );
  } else {
    return 1 + depth - minDepth;
  }
}

//
// TransformSVG
// 

real
TransformSVG::rounded( real x ) const
{
  return round( 100*x ) / 100.0f;
} 

real
TransformSVG::mapY( real y ) const
{
  const real ppmm = 720/254.0f;
  return rounded( ( 297 * ppmm ) - ( y * _scale + _deltaY ) );
}

real
TransformSVG::mapWidth( real width ) const
{
  const real ppmm = 72.0/25.4f;
  return round( 1000 * width / ppmm ) / 1000.0;
}

void
TransformSVG::setBoundingBox( const Rect & rect )
{
  const real ppmm = 720/254.0f;
  if ( rect.height/rect.width > 27.7/19.0 ) {
    _scale = 277 * ppmm / rect.height;
  } else {
    _scale = 190 * ppmm / rect.width;
  }
  _deltaX = 0.5 * 210 * ppmm - _scale * ( rect.left + 0.5 * rect.width );
  _deltaY = 0.5 * 297 * ppmm - _scale * ( rect.top - 0.5 * rect.height );
}
