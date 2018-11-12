#include "rect.h"
#include "shapes.h"
#include <cmath>
#include <vector>
#include <sstream>
#include <cstring>

bool shapeGreaterDepth( const Form *s1, const Form *s2 )
{
  return s1->depth > s2->depth;
}

std::string
Form::svgProperties( const TransformSVG & transform ) const
{
  static const char * capStrings[3] = { "butt", "round", "square" };
  static const char * joinStrings[3] = { "miter", "round", "bevel" };
  std::stringstream str;
  if ( penColor != Color::none ) {
    str << " fill=\"" << fillColor.svg() << '"'
	<< " stroke=\"" << penColor.svg() << '"'
	<< " stroke-width=\"" << transform.mapWidth( lineWidth ) << "mm\""
	<< " style=\"stroke-linecap:" << capStrings[ lineCap ]  
	<< ";stroke-linejoin:" << joinStrings[ lineJoin ] << '"'
	<< fillColor.svgAlpha( " fill" )
	<< penColor.svgAlpha( " stroke" );
  } else  {
    str << " fill=\"" << fillColor.svg() << '"'
	<< " stroke=\"" << fillColor.svg() << '"'
	<< " stroke-width=\"0.5px\""
	<< " style=\"stroke-linecap:round;stroke-linejoin:round;\""
	<< fillColor.svgAlpha( " fill" )
	<< fillColor.svgAlpha( " stroke" );
  }
  return str.str();
}

void
Line::flushPostscript( std::ostream & stream,
		       const TransformEPS & transform ) const
{
  stream << "n "
	 << transform.mapX( x1 ) << " " 
	 << transform.mapY( y1 ) << " " 
	 << "m "
	 << transform.mapX( x2 ) << " " 
	 << transform.mapY( y2 ) << " " 
	 << "l " << penColor.postscript() << " srgb stroke" << std::endl;
}

void
Line::flushFIG( std::ostream & stream,
		const TransformFIG & transform,
		std::map<Color,int> & colormap ) const
{
  stream << "2 1 0 ";
  // Thickness
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color
  stream << colormap[ penColor ] << " ";
  // Fill color
  stream << "0 ";
  // Depth
  stream << transform.mapDepth( depth ) << " ";
  // Pen style
  stream <<  "-1 ";
  // Area fill, style val, join style, cap style, radius, f_arrow, b_arrow
  stream << "-1 0.000 " << lineJoin << " " << lineCap << " -1 0 0 ";
  // Number of points
  stream << "2\n";
  stream << "         ";
  stream << static_cast<int>( transform.mapX( x1 ) ) << " " 
	 << static_cast<int>( transform.mapY( y1 ) ) << " "
	 << static_cast<int>( transform.mapX( x2 ) ) << " " 
	 << static_cast<int>( transform.mapY( y2 ) ) << std::endl;
}

void
Line::flushSVG( std::ostream & stream,
		const TransformSVG & transform ) const
{
  stream << "<line x1=\"" << transform.mapX( x1 ) << "\""
	 << " y1=\"" << transform.mapY( y1 ) << "\""
	 << " x2=\"" << transform.mapX( x2 ) << "\""
	 << " y2=\"" << transform.mapY( y2 ) << "\""
    	 << svgProperties( transform ) 
	 << " />" << std::endl;
}

Rect
Line::boundingBox() const
{
  Rect rect;
  if ( x1 > x2 ) {
    rect.width = x1 - x2;
    rect.left = x2;
  } else {
    rect.width = x2 - x1;
    rect.left = x1;
  }
  if ( y1 > y2 ) {
    rect.top = y1;
    rect.height = y1 - y2;
  } else {
    rect.top = y2;
    rect.height = y2 - y1;
  }
  return rect;
}

void
Arrow::flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const
{
  real dx = x1 - x2;
  real dy = y1 - y2;
  real norm = sqrt( dx*dx + dy*dy );
  dx /= norm;
  dy /= norm;
  dx *= 10*lineWidth;
  dy *= 10*lineWidth;
  real ndx1 = dx*cos(0.3)-dy*sin(0.3);
  real ndy1 = dx*sin(0.3)+dy*cos(0.3);
  real ndx2 = dx*cos(-0.3)-dy*sin(-0.3);
  real ndy2 = dx*sin(-0.3)+dy*cos(-0.3);
    
  stream << "n "
	 << transform.mapX( x1 ) << " " 
	 << transform.mapY( y1 ) << " " 
	 << "m "
	 << transform.mapX( x2 ) << " " 
	 << transform.mapY( y2 ) << " " 
	 << "l stroke" << std::endl;

  if ( filled() ) { 
    stream << "n "
	   << transform.mapX( x2 ) + ndx1 << " " 
	   << transform.mapY( y2 ) + ndy1 << " " 
	   << "m "
	   << transform.mapX( x2 ) << " " 
	   << transform.mapY( y2 ) << " l " 
	   << transform.mapX( x2 ) + ndx2 << " " 
	   << transform.mapY( y2 ) + ndy2 << " ";
    stream  << "l cp " << penColor.postscript() << " srgb  fill" << std::endl;
  }
  
  stream << "n "
	 << transform.mapX( x2 ) + ndx1 << " " 
	 << transform.mapY( y2 ) + ndy1 << " " 
	 << "m "
	 << transform.mapX( x2 ) << " " 
	 << transform.mapY( y2 ) << " l " 
	 << transform.mapX( x2 ) + ndx2 << " " 
	 << transform.mapY( y2 ) + ndy2 << " l"
	 << " " << penColor.postscript() << " srgb stroke" << std::endl;
}

void
Arrow::flushFIG( std::ostream & stream,
		const TransformFIG & transform,
		std::map<Color,int> & colormap ) const
{
  stream << "2 1 0 ";
  // Thickness
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color
  stream << colormap[ penColor ] << " ";
  // Fill color
  stream << colormap[ penColor ] << " ";
  // Depth
  stream << transform.mapDepth( depth ) << " ";
  // Pen style
  stream <<  "-1 ";
  // Area fill, style val, join style, cap style, radius, f_arrow, b_arrow
  stream << "-1 0.000 " << lineJoin << " " << lineCap << " -1 1 0 ";
  // Number of points
  stream << "2\n";
  if ( filled() ) 
    stream << "       1 1 1.00 60.00 120.00\n";
  else 
    stream << "       0 0 1.00 60.00 120.00\n";
  stream << "         ";
  stream << static_cast<int>( transform.mapX( x1 ) ) << " " 
	 << static_cast<int>( transform.mapY( y1 ) ) << " "
	 << static_cast<int>( transform.mapX( x2 ) ) << " " 
	 << static_cast<int>( transform.mapY( y2 ) ) << std::endl;
}

void
Arrow::flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const
{
  real dx = x1 - x2;
  real dy = y1 - y2;
  real norm = sqrt( dx*dx + dy*dy );
  dx /= norm;
  dy /= norm;
  dx *= 10*lineWidth;
  dy *= 10*lineWidth;
  real ndx1 = dx*cos(0.3)-dy*sin(0.3);
  real ndy1 = dx*sin(0.3)+dy*cos(0.3);
  real ndx2 = dx*cos(-0.3)-dy*sin(-0.3);
  real ndy2 = dx*sin(-0.3)+dy*cos(-0.3);
  
  stream << "<g>" << std::endl;
  // The line
  stream << " <path "
	 << "d=\"M " << transform.mapX(x1) << " " << transform.mapY(y1)
	 << " L " << transform.mapX(x2) << " " << transform.mapY(y2) << " z\""
	 << " fill=\"none\" stroke=\"" << penColor.svg() << "\""
    	 << penColor.svgAlpha( " stroke" )
	 << " stroke-width=\"" << transform.mapWidth( lineWidth ) << "mm\" />";
  
  // The arrow
  stream << " <polygon";
  stream << " fill=\"" << fillColor.svg() << "\"";
  stream << " stroke=\"" << penColor.svg() << "\""
	 << " stroke-width=\"" << transform.mapWidth( 0.33 * lineWidth ) << "mm\""
	 << " style=\"stroke-linecap:butt;stroke-linejoin:miter\""
	 << fillColor.svgAlpha( " fill" )
	 << penColor.svgAlpha( " stroke" )
	 << " points=\""
	 << transform.mapX( x2 ) + ndx1 << "," 
	 << transform.mapY( y2 ) - ndy1 << " "
	 << transform.mapX( x2 ) << "," 
	 << transform.mapY( y2 ) << " " 
	 << transform.mapX( x2 ) + ndx2 << "," 
	 << transform.mapY( y2 ) - ndy2 << " "
	 << transform.mapX( x2 ) + ndx1 << "," 
	 << transform.mapY( y2 ) - ndy1 << "\" />" << std::endl;
  stream << "</g>" << std::endl;
}

void
Rectangle::flushPostscript( std::ostream & stream,
			    const TransformEPS & transform ) const
{
  if ( filled() ) {
    stream << "n " 
	   << transform.mapX( x ) << " "  << transform.mapY( y ) << " m " 
	   << transform.mapX( x + width ) << " " << transform.mapY( y ) << " l "
	   << transform.mapX( x + width ) << " " << transform.mapY( y - height ) << " l "
	   << transform.mapX( x ) << " " << transform.mapY( y - height ) << " l cp";
    stream << " " << fillColor.postscript() << " srgb";
    stream << " fill" << std::endl;  
  }  
  if ( penColor != Color::none ) {
    stream << "n " 
	   << transform.mapX( x ) << " "  << transform.mapY( y ) << " m " 
	   << transform.mapX( x + width ) << " " << transform.mapY( y ) << " l "
	   << transform.mapX( x + width ) << " " << transform.mapY( y - height ) << " l "
	   << transform.mapX( x ) << " " << transform.mapY( y - height ) << " l cp";
    stream << " " << penColor.postscript() << " srgb";
    stream << " stroke" << std::endl;
  }
}

void
Rectangle::flushFIG( std::ostream & stream,
		     const TransformFIG & transform,
		     std::map<Color,int> & colormap ) const
{
  stream << "2 2 0 ";
  // Thickness
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color
  stream << colormap[ penColor ] << " ";
  // Fill color
  stream << colormap[ fillColor ] << " ";
  // Depth
  stream << transform.mapDepth( depth ) << " ";
  // Pen style
  stream <<  "-1 ";
  // Area fill, style val, join style, cap style, radius, f_arrow, b_arrow, number of points
  if ( filled() ) 
    stream << "20 0.000 " << lineJoin << " " << lineCap << " -1 0 0 5\n";
  else
    stream << "-1 0.000 " << lineJoin << " " << lineCap << " -1 0 0 5\n";
  stream << "         ";
  
  stream << static_cast<int>( transform.mapX( x ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << " "
	 << static_cast<int>( transform.mapX( x + width ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << " "
	 << static_cast<int>( transform.mapX( x + width ) ) << " "
	 << static_cast<int>( transform.mapY( y - height ) ) << " "
	 << static_cast<int>( transform.mapX( x ) ) << " "
	 << static_cast<int>( transform.mapY( y - height ) ) << " "
	 << static_cast<int>( transform.mapX( x ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << "\n";
}

void
Rectangle::flushSVG( std::ostream & stream,
		     const TransformSVG & transform ) const
{

  stream << "<rect x=\"" << transform.mapX( x ) << '"'
	 << " y=\"" << transform.mapY( y )  << '"'
	 << " width=\"" << transform.scale( width ) << '"'
	 << " height=\"" << transform.scale( height ) << '"'
	 << svgProperties( transform ) 
	 << " />" << std::endl;
}

Rect
Rectangle::boundingBox() const
{
  Rect rect;
  if ( width > 0.0f ) { 
    rect.left = x;
    rect.width = width;
  } else {
    rect.left = x + width;
    rect.width = fabs( width );
  }
  if ( height > 0.0f ) { 
    rect.top = y;
    rect.height = height;
  } else {
    rect.top = y + height;
    rect.height = fabs( height );
  }
  return rect;
}

void
Circle::flushPostscript( std::ostream & stream,
			 const TransformEPS & transform ) const
{
  if ( filled() ) {
    stream << "n " << transform.mapX( x + radius ) << " "  << transform.mapY( y ) << " m " 
	   << transform.mapX( x ) << " " << transform.mapY( y ) << " "
	   << transform.scale( radius ) << " 0.0 360.0 arc ";
    stream << " " << fillColor.postscript() << " srgb";
    stream << " fill " << std::endl;  
  }
  
  if ( penColor != Color::none ) {
    stream << "n " << transform.mapX( x + radius ) << " "  << transform.mapY( y ) << " m " 
	   << transform.mapX( x ) << " " << transform.mapY( y ) << " "
	   << transform.scale( radius ) << " 0.0 360.0 arc ";
    stream << " " << penColor.postscript() << " srgb";
    stream << " stroke " << std::endl;
  }
}

void
Circle::flushFIG( std::ostream & stream,
		  const TransformFIG & transform,
		  std::map<Color,int> & colormap ) const
{
  // Ellipse, Sub type, Line style, Thickness
  stream << "1 3 0 ";
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color, Fill color
  stream << colormap[ penColor ] << " " << colormap[ fillColor ] << " ";
  // Depth, Pen style, Area fill, Style val, Direction 
  if ( filled() )
    stream << transform.mapDepth( depth ) << " -1 20 0.000 1 0.000 ";
  else
    stream << transform.mapDepth( depth ) << " -1 -1 0.000 1 0.000 ";
  stream << static_cast<int>( transform.mapX( x ) ) << " " 
	 << static_cast<int>( transform.mapY( y ) ) << " " 
	 << static_cast<int>( transform.scale( radius ) ) << " " 
	 << static_cast<int>( transform.scale( radius ) ) << " " 
	 << static_cast<int>( transform.mapX( x ) ) << " " 
	 << static_cast<int>( transform.mapY( y ) ) << " " 
	 << static_cast<int>( transform.mapX( x ) + transform.scale( radius ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << "\n";
}

void
Circle::flushSVG( std::ostream & stream,
		  const TransformSVG & transform ) const
{
  stream << "<circle cx=\"" << transform.mapX( x ) << '"'
	 << " cy=\"" << transform.mapY( y ) << '"'
	 << " r=\"" << transform.scale( radius ) << '"'
	 << svgProperties( transform ) 
	 << " />" << std::endl;
}

Rect
Circle::boundingBox() const
{
  return Rect( x - radius, y + radius, 2 * radius, 2 * radius ); 
}

void
Ellipse::flushPostscript( std::ostream & stream,
			 const TransformEPS & transform ) const
{
  real yScale = yRadius / xRadius;
  
  if ( filled() ) {
    stream << "gs " << transform.mapX( x ) << " " << transform.mapY( y ) << " tr";
    stream << " " << 1.0 << " " << yScale << " sc";
    stream << " n " << transform.scale( xRadius ) << " 0 m " 
	   << " 0 0 " << transform.scale( xRadius ) << " 0.0 360.0 arc ";
    stream << " " << fillColor.postscript() << " srgb";
    stream << " fill gr" << std::endl;  
  }
  
  if ( penColor != Color::none ) {
    stream << "gs " << transform.mapX( x ) << " " << transform.mapY( y ) << " tr";
    stream << " " << 1.0 << " " << yScale << " sc";
    stream << " n " << transform.scale( xRadius ) << " 0 m " 
	   << " 0 0 " << transform.scale( xRadius ) << " 0.0 360.0 arc ";
    stream << " " << penColor.postscript() << " srgb";
    stream << " stroke gr" << std::endl;  
  }
}

void
Ellipse::flushFIG( std::ostream & stream,
		  const TransformFIG & transform,
		  std::map<Color,int> & colormap ) const
{
  // Ellipse, Sub type, Line style, Thickness
  stream << "1 3 0 ";
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color, Fill color
  stream << colormap[ penColor ] << " " << colormap[ fillColor ] << " ";
  // Depth, Pen style, Area fill, Style val, Direction 
  if ( filled() )
    stream << transform.mapDepth( depth ) << " -1 20 0.000 1 0.000 ";
  else
    stream << transform.mapDepth( depth ) << " -1 -1 0.000 1 0.000 ";
  stream << static_cast<int>( transform.mapX( x ) ) << " " 
	 << static_cast<int>( transform.mapY( y ) ) << " " 
	 << static_cast<int>( transform.scale( xRadius ) ) << " " 
	 << static_cast<int>( transform.scale( yRadius ) ) << " " 
	 << static_cast<int>( transform.mapX( x ) ) << " " 
	 << static_cast<int>( transform.mapY( y ) ) << " " 
	 << static_cast<int>( transform.mapX( x ) + transform.scale( xRadius ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << "\n";
}

void
Ellipse::flushSVG( std::ostream & stream,
		  const TransformSVG & transform ) const
{
  stream << "<ellipse cx=\"" << transform.mapX( x ) << '"'
	 << " cy=\"" << transform.mapY( y ) << '"'
	 << " rx=\"" << transform.scale( xRadius ) << '"'
	 << " ry=\"" << transform.scale( yRadius ) << '"'
    	 << svgProperties( transform ) 
	 << " />" << std::endl;
}

Rect
Ellipse::boundingBox() const
{
  return Rect( x - xRadius, y + yRadius, 2 * xRadius, 2 * yRadius ); 
}

void
Polyline::flushPostscript( std::ostream & stream,
			   const TransformEPS & transform ) const
{
  std::vector<Point>::const_iterator i = points.begin();
  std::vector<Point>::const_iterator end = points.end();
  
  if ( filled() ) {
    stream << "n " << transform.mapX( i->x ) << " " << transform.mapY( i->y ) << " m";
    ++i;
    while ( i != end ) {
      stream << " " << transform.mapX( i->x ) << " " << transform.mapY( i->y ) << " l";
      ++i;
    }
    if ( closed ) stream << " cp";
    stream << " " << fillColor.postscript() << " srgb";
    stream << " fill" << std::endl;  
  }

  i = points.begin();
  if ( penColor != Color::none ) {
    stream << "n " << transform.mapX( i->x ) << " " << transform.mapY( i->y ) << " m";
    ++i;
    while ( i != end ) {
      stream << " " << transform.mapX( i->x ) << " " << transform.mapY( i->y ) << " l";
      ++i;
    }
    if ( closed ) stream << " cp";
    stream << " " << penColor.postscript() << " srgb";
    stream << " stroke" << std::endl;
  }
}

void
Polyline::flushFIG( std::ostream & stream,
		    const TransformFIG & transform,
		    std::map<Color,int> & colormap ) const
{
  if ( closed ) 
    stream << "2 3 0 ";
  else
    stream << "2 1 0 ";
  // Thickness
  stream << ( penColor.valid()?transform.mapWidth( lineWidth ):0 ) << " ";
  // Pen color
  stream << colormap[ penColor ] << " ";
  // Fill color
  stream << colormap[ fillColor ] << " ";
  // Depth
  stream << transform.mapDepth( depth ) << " ";
  // Pen style
  stream <<  "-1 ";
  // Area fill, style val, join style, cap style, radius, f_arrow, b_arrow
  if ( filled() ) 
    stream << "20 0.000 " << lineJoin << " " << lineCap << " -1 0 0 ";
  else
    stream << "-1 0.000 " << lineJoin << " " << lineCap << " -1 0 0 ";
  // Number of points
  stream << points.size() + closed << std::endl;
  stream << "         ";

  std::vector<Point>::const_iterator i = points.begin();
  std::vector<Point>::const_iterator end = points.end();
  while ( i != end ) {
    stream << " " << static_cast<int>( transform.mapX( i->x ) )
	   << " " << static_cast<int>( transform.mapY( i->y ) );
    ++i;
  }
  if ( closed ) { 
    stream << " " << static_cast<int>( transform.mapX( points.begin()->x ) )
	   << " " << static_cast<int>( transform.mapY( points.begin()->y ) );
  }  
  stream << std::endl;
}

void
Polyline::flushSVG( std::ostream & stream,
		    const TransformSVG & transform ) const
{
  std::vector<Point>::const_iterator i = points.begin();
  std::vector<Point>::const_iterator end = points.end();
  int count = 0;
  if ( closed )
    stream << "<polygon ";
  else
    stream << "<polyline";
  
  stream << svgProperties( transform ) << std::endl;
  stream << "          points=\"";

  stream << transform.mapX( i->x ) << "," << transform.mapY( i->y );
  ++i;
  while ( i != end ) {
    stream << " " << transform.mapX( i->x ) << "," << transform.mapY( i->y );
    ++i;
    count = ( count + 1 ) % 6;
    if ( !count ) stream << "\n                  ";
  }
  stream << "\" />" << std::endl;
}

Rect
Polyline::boundingBox() const
{
  Rect rect;
  std::vector< Point >::const_iterator i = points.begin();
  std::vector< Point >::const_iterator end = points.end();
  rect.top = i->y;
  rect.left = i->x;
  rect.width = 0.0;
  rect.height = 0.0;
  ++i;
  while ( i != end ) { 
    if ( i->x < rect.left ) { 
      real dw = rect.left - i->x;
      rect.left = i->x;
      rect.width += dw;
    } else if ( i->x > rect.left + rect.width ) {
      rect.width = i->x - rect.left;
    }
    if ( i->y > rect.top ) { 
      real dh = i->y - rect.top;
      rect.top = i->y;
      rect.height += dh;
    } else if ( i->y < rect.top - rect.height ) {
      rect.height = rect.top - i->y;
    }
    ++i;
  }
  return rect;
}


GouraudTriangle::GouraudTriangle( const Point & p0, const Color & color0,
				  const Point & p1, const Color & color1,
				  const Point & p2, const Color & color2,
				  int subdivisions,
				  unsigned int depth )
  : Polyline( std::vector<Point>(), true, Color::none, Color::none,
	      0.0f, ButtCap, MiterJoin, depth ),
    color0( color0 ), color1( color1 ), color2( color2 ), subdivisions( subdivisions ) {
  points.push_back( p0 );
  points.push_back( p1 );
  points.push_back( p2 );

  Form::fillColor.red( ( color0.red() + color1.red() + color2.red() ) / 3 );
  Form::fillColor.green( ( color0.green() + color1.green() + color2.green() ) / 3 );
  Form::fillColor.blue( ( color0.blue() + color1.blue() + color2.blue() ) / 3 );
}

GouraudTriangle::GouraudTriangle( const Point & p0, real brightness0,
				  const Point & p1, real brightness1,
				  const Point & p2, real brightness2,
				  const Color & fillColor,
				  int subdivisions,
				  unsigned int depth )
  : Polyline( std::vector<Point>(), true, Color::none, Color::none,
	      0.0f, ButtCap, MiterJoin, depth ),
    color0( fillColor ), color1( fillColor ), color2( fillColor ), subdivisions( subdivisions )
{
  points.push_back( p0 );
  points.push_back( p1 );
  points.push_back( p2 );
  color0.red( static_cast<unsigned char>( std::min( 255.0, color0.red() * brightness0 ) ) );
  color0.green( static_cast<unsigned char>( std::min( 255.0, color0.green() * brightness0 ) ) );
  color0.blue( static_cast<unsigned char>( std::min( 255.0, color0.blue() * brightness0 ) ) );
  color1.red( static_cast<unsigned char>( std::min( 255.0, color1.red() * brightness1 ) ) );
  color1.green( static_cast<unsigned char>( std::min( 255.0, color1.green() * brightness1 ) ) );
  color1.blue( static_cast<unsigned char>( std::min( 255.0, color1.blue() * brightness1 ) ) );
  color2.red( static_cast<unsigned char>( std::min( 255.0, color2.red() * brightness2 ) ) );
  color2.green( static_cast<unsigned char>( std::min( 255.0, color2.green() * brightness2 ) ) );
  color2.blue( static_cast<unsigned char>( std::min( 255.0, color2.blue() * brightness2 ) ) );

  Form::fillColor.red( ( color0.red() + color1.red() + color2.red() ) / 3 );
  Form::fillColor.green( ( color0.green() + color1.green() + color2.green() ) / 3 );
  Form::fillColor.blue( ( color0.blue() + color1.blue() + color2.blue() ) / 3 );
}
 
void
GouraudTriangle::flushPostscript( std::ostream & stream,
				  const TransformEPS & transform ) const
{
  if ( ! subdivisions ) {
    Polyline::flushPostscript( stream, transform );
    return;
  }
  const Point & p0 = points[0];
  const Point & p1 = points[1];
  const Point & p2 = points[2];  
  Point p01( 0.5*(p0.x+p1.x), 0.5*(p0.y+p1.y) );
  Color c01( (color0.red() + color1.red())/2, 
	     (color0.green() + color1.green())/2, 
	     (color0.blue() + color1.blue())/2 );
  Point p12( 0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y) );
  Color c12( (color1.red() + color2.red())/2, 
	     (color1.green() + color2.green())/2, 
	     (color1.blue() + color2.blue())/2 );
  Point p20( 0.5*(p2.x+p0.x), 0.5*(p2.y+p0.y) );
  Color c20( (color2.red() + color0.red())/2, 
	     (color2.green() + color0.green())/2, 
	     (color2.blue() + color0.blue())/2 );
  GouraudTriangle( p0, color0, p20, c20, p01, c01, subdivisions - 1, depth ).flushPostscript( stream, transform );
  GouraudTriangle( p1, color1, p01, c01, p12, c12, subdivisions - 1, depth ).flushPostscript( stream, transform );
  GouraudTriangle( p2, color2, p20, c20, p12, c12, subdivisions - 1, depth ).flushPostscript( stream, transform );
  GouraudTriangle( p01, c01, p12, c12, p20, c20,  subdivisions - 1, depth ).flushPostscript( stream, transform );
}


void
GouraudTriangle::flushFIG( std::ostream & stream,
			   const TransformFIG & transform,
			   std::map<Color,int> & colormap ) const
{
  Polyline::flushFIG( stream, transform, colormap );
}

void
GouraudTriangle::flushSVG( std::ostream & stream,
			   const TransformSVG & transform ) const
{
  if ( ! subdivisions ) {
    Polyline::flushSVG( stream, transform );
    return;
  }
  const Point & p0 = points[0];
  const Point & p1 = points[1];
  const Point & p2 = points[2];  
  Point p01( 0.5*(p0.x+p1.x), 0.5*(p0.y+p1.y) );
  Color c01( (color0.red() + color1.red())/2, 
	     (color0.green() + color1.green())/2, 
	     (color0.blue() + color1.blue())/2 );
  Point p12( 0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y) );
  Color c12( (color1.red() + color2.red())/2, 
	     (color1.green() + color2.green())/2, 
	     (color1.blue() + color2.blue())/2 );
  Point p20( 0.5*(p2.x+p0.x), 0.5*(p2.y+p0.y) );
  Color c20( (color2.red() + color0.red())/2, 
	     (color2.green() + color0.green())/2, 
	     (color2.blue() + color0.blue())/2 );
  GouraudTriangle( p0, color0, p20, c20, p01, c01, subdivisions - 1, depth ).flushSVG( stream, transform );
  GouraudTriangle( p1, color1, p01, c01, p12, c12, subdivisions - 1, depth ).flushSVG( stream, transform );
  GouraudTriangle( p2, color2, p20, c20, p12, c12, subdivisions - 1, depth ).flushSVG( stream, transform );
  GouraudTriangle( p01, c01, p12, c12, p20, c20,  subdivisions - 1, depth ).flushSVG( stream, transform ); 
}

void
Text::flushPostscript( std::ostream & stream,
		       const TransformEPS & transform ) const
{
  stream << "/" << font << " ff " << size << " scf sf";
  stream << " " << transform.mapX( x ) << " " << transform.mapY( y ) << " m";
  stream << " (" << text << ")" 
	 << " " << penColor.postscript() << " srgb"
	 << " sh" << std::endl;
}

void
Text::flushFIG( std::ostream & stream,
		const TransformFIG & transform,
		std::map<Color,int> & colormap ) const
{
  extern const char * XFigPostscriptFontnames[];
  int i = 0;
  for ( i = 0 ; XFigPostscriptFontnames[ i ] && strcmp( XFigPostscriptFontnames[ i ], 
							font.c_str() ); ++i );
  if ( ! XFigPostscriptFontnames[ i ] ) i = 0;
  stream << "4 0 " ;
  // Color, depth, unused, Font
  stream << colormap[ penColor ] <<  " " << transform.mapDepth( depth ) << " -1 " << i << " ";
  // Font size, Angle, Font flags
  stream << size << " 0.000 4 ";
  // Height
  stream << static_cast<int>( size * 135 / 12.0 ) << " ";
  // Width
  stream << static_cast<int>( text.size() * size * 135 / 12.0 ) << " ";
  // x y 
  stream << static_cast<int>( transform.mapX( x ) ) << " "
	 << static_cast<int>( transform.mapY( y ) ) << " ";
  stream << text << "\\001\n";
}

void
Text::flushSVG( std::ostream & stream,
		       const TransformSVG & transform ) const
{
  stream << "<text x=\"" << transform.mapX( x )
	 << "\" y=\"" << transform.mapY( y ) << "\" "
	 << " font-family=\"" << font << "\""
	 << " font-size=\"" << size << "\""
	 << " fill=\"" << penColor.svg() 
	 << fillColor.svgAlpha( " fill" )
	 << penColor.svgAlpha( " stroke" )
	 << "\" >" << std::endl
	 << text << std::endl
	 << "</text>" << std::endl;
}

Rect
Text::boundingBox() const
{
  return Rect( x, y, 0, 0 );
}
