#ifndef _BOARD_COLOR_H_
#define _BOARD_COLOR_H_

#include <ostream>
#include <string>

#include "../Global/global_precision.h"

/**
 * Color structure.
 * @brief Structure representing an RGB triple.
 */
class Color {

public:
  Color( unsigned char red, unsigned char  green, unsigned char  blue, unsigned char alpha = 255 )
    :_red(red),_green(green),_blue(blue),_alpha(alpha) { }

  Color( unsigned char gray, unsigned char alpha = 255 )
    :_red(gray),_green(gray),_blue(gray),_alpha(alpha) { }

  Color( const bool valid = true )
    :_red(-1),_green(-1),_blue(-1),_alpha(255)
  { 
    if ( valid ) {
      _red = _green = _blue = 0;
    }
  }
    
  inline void red( unsigned char red );
  inline void green( unsigned char green );
  inline void blue( unsigned char blue );
  inline void alpha( unsigned char alpha );

  inline unsigned char red() const;
  inline unsigned char green() const;
  inline unsigned char blue() const;
  inline unsigned char alpha() const;
  
  inline Color & setRGBi( const unsigned char red,
			  const unsigned char green,
			  const unsigned char blue,
			  const unsigned char alpha = 255 );
  
  Color & setRGBf( real red, 
		   real green,
		   real blue,
		   real alpha = 1.0 );
  
  bool operator==( const Color & other ) const;

  bool operator!=( const Color & other ) const;

  bool operator<( const Color & other ) const;

  void flushPostscript( std::ostream & ) const;

  std::string svg() const;

  /** 
   * Return a an SVG parameter string for the opacity value.
   * 
   * @param prefix A prefix string to be appended to the returned 
   * string if not empty.
   * 
   * @return An empty string if alpha == 255, otherwise the string <prefix>-opacity="<alpha-value>".
   */
  std::string svgAlpha( const char * prefix ) const;

  std::string postscript() const;
  
  inline bool valid() const { return (*this) != Color::none; }

public:
  static const Color none;
  static const Color black;
  static const Color white;

private:
  int _red;			/**< The red component. */
  int _green;			/**< The green component. */
  int _blue;			/**< The blue component. */
  int _alpha;			/**< The opacity. */
};

inline Color &
Color::setRGBi( const unsigned char red,
		const unsigned char green,
		const unsigned char blue,
		const unsigned char alpha ) {
  _red = red;
  _green = green;
  _blue = blue;
  _alpha = alpha;
  return *this;
}


inline void
Color::red( const unsigned char red )
{
  _red = red;
}

inline void
Color::green( unsigned char green )
{
  _green = green;
}

inline void 
Color::blue( unsigned char blue )
{
  _blue = blue;
}

inline void 
Color::alpha( unsigned char alpha )
{
  _alpha = alpha;
}

inline unsigned char Color::red() const { return _red; } 
inline unsigned char Color::green() const { return _green; } 
inline unsigned char Color::blue() const { return _blue; } 
inline unsigned char Color::alpha() const { return _alpha; } 

#endif // _COLOR_H_

/*-----------------------------------------------------------------------------+
 '$Id: color.h,v 1.3 2008/11/17 19:23:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
