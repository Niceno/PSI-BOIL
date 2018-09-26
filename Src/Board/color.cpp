#include "color.h"
#include <cstdio>
using std::string;

const Color Color::none(false);
const Color Color::black((unsigned char)0,(unsigned char)0,(unsigned char)0);
const Color Color::white((unsigned char)255,(unsigned char)255,(unsigned char)255);

Color &
Color::setRGBf( real red,
		real green,
		real blue,
		real alpha  ) {
  if ( red > 1.0f ) red = 1.0f;
  if ( red < 0.0f ) red = 0.0f;
  _red = static_cast<unsigned char>( 255 * red );
  if ( green > 1.0f ) green = 1.0f;
  if ( green < 0.0f ) green = 0.0f;
  _green = static_cast<unsigned char>( 255 * green );
  if ( blue > 1.0f ) blue = 1.0f;
  if ( blue < 0.0f ) blue = 0.0f;
  _blue = static_cast<unsigned char>( 255 * blue );
  if ( alpha > 1.0f ) alpha = 1.0f;
  if ( alpha < 0.0f ) alpha = 0.0f;
  _alpha = static_cast<unsigned char>( 255 * alpha );
  return *this;
}

bool
Color::operator==( const Color & other ) const
{
  return _red == other._red  && _green == other._green && _blue == other._blue
    && _alpha == other._alpha;
}

bool
Color::operator!=( const Color & other ) const
{
  return _red != other._red  || _green != other._green || _blue != other._blue
    || _alpha != other._alpha;
}

bool
Color::operator<( const Color & other ) const
{
  if ( _red < other._red ) return true;
  if ( _red == other._red ) {
    if ( _green < other._green ) return true;
    if ( _green == other._green ) { 
      if ( _blue < other._blue ) return true;
      if ( _blue == other._blue )
	return _alpha < other._alpha;
    }
  }
  return false;
}

void
Color::flushPostscript( std::ostream & stream ) const
{
  stream << (_red/255.0) << " "
	 << (_green/255.0) << " "
	 << (_blue/255.0) << " srgb\n";
}

string
Color::postscript() const
{
  char buffer[255];
  sprintf( buffer, "%.4f %.4f %.4f", _red/255.0, _green/255.0, _blue/255.0 );
  return buffer;
}

string
Color::svg() const
{
  char buffer[255];
  if ( *this == Color::none ) return "none";
  sprintf( buffer, "rgb(%d,%d,%d)", _red, _green, _blue );
  return buffer;
}

string
Color::svgAlpha( const char * prefix ) const
{
  char buffer[255];
  if ( _alpha == 255 || *this == Color::none ) return "";
  sprintf( buffer, " %s-opacity=\"%f\"", prefix, _alpha/255.0f );
  return buffer;
}

/*-----------------------------------------------------------------------------+
 '$Id: color.cpp,v 1.3 2010/01/20 10:19:08 sato Exp $'/
+-----------------------------------------------------------------------------*/
