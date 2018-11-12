#ifndef _BOARD_BOARD_H_
#define _BOARD_BOARD_H_

#include <ostream>
#include <string>

#include <vector>

#include "../Global/global_precision.h"
#include "point.h"
#include "shapes.h"

/** 
 * The Board class.
 * @brief Class for EPS, FIG or SVG drawings.
 * @version 0.5
 */
class Board {
 public:

  /** 
   * Constructs a new board and sets the background color, if any.
   * 
   * @param backgroundColor A color for the drawing's background.
   */
  Board( const Color & backgroundColor = Color::none );

  ~Board();
  
  /** 
   * Clears the board with a given background color.
   * 
   * @param color The board background color (may be Color::none).
   */
  void clear( const Color & color = Color::none );

  /** 
   * Clears the board and set the background color from an RGB triple.
   * 
   * @param red 
   * @param green 
   * @param blue 
   */
  inline void clear( unsigned char red, unsigned char green, unsigned char blue );

  /** 
   * Draws a line from (x1,y1) to (x2,y2).
   * 
   * @param x1 First coordinate of the first extremity.
   * @param y1 Second coordinate of the first extremity.
   * @param x2 First coordinate of the second extremity.
   * @param y2 Second coordinate of the second extremity.
   * @param depth Depth of the line.
   */
  void drawLine( real x1, real y1, real x2, real y2, 
		 unsigned int depth = 0 );

  /** 
   * Draws a line from (x1,y1) to (x2,y2) with an arrow at (x2,y2).
   * 
   * @param x1 First coordinate of the first extremity.
   * @param y1 Second coordinate of the first extremity.
   * @param x2 First coordinate of the second extremity.
   * @param y2 Second coordinate of the second extremity.
   * @param filled Whether or not the arrow is filled.
   * @param depth Depth of the line.
   */
  void drawArrow( real x1, real y1, real x2, real y2,
		  bool filled = false,
		  unsigned int depth = 0 );

  /** 
   * Draws a triangle.
   * 
   * @param x1 First coordinate of the first vertex.
   * @param y1 Second coordinate of the first vertex.
   * @param x2 First coordinate of the second vertex.
   * @param y3 Second coordinate of the second vertex.
   * @param x3 First coordinate of the third vertex.
   * @param y3 Second coordinate of the third vertex.
   * @param depth Depth of the triangle.
   */
  void drawTriangle( real x1, real y1, 
		     real x2, real y2, 
		     real x3, real y3, 
		     unsigned int depth = 0 );

  /** 
   * Draws a triangle.
   * 
   * @param p1 First vertex.
   * @param p2 Second vertex.
   * @param p3 Third vertex.
   * @param depth Depth of the triangle.
   */
  void drawTriangle( const Point & p1,
		     const Point & p2, 
		     const Point & p3, 
		     unsigned int depth = 0 );

  /** 
   * Draws a filled triangle.
   * 
   * @param x1 First coordinate of the first vertex.
   * @param y1 Second coordinate of the first vertex.
   * @param x2 First coordinate of the second vertex.
   * @param y3 Second coordinate of the second vertex.
   * @param x3 First coordinate of the third vertex.
   * @param y3 Second coordinate of the third vertex.
   * @param depth Depth of the triangle.
   */
  void fillTriangle( real x1, real y1, 
		     real x2, real y2, 
		     real x3, real y3, 
		     unsigned int depth = 0 );

  /** 
   * Draws a triangle with Gouraud-like shaded colors. 
   * 
   * @param p1 
   * @param color1 
   * @param p2 
   * @param color2 
   * @param p3 
   * @param color3 
   * @param divisions number of triangle subdivisions.
   * @param depth The depth of the triangle.
   */
  void fillGouraudTriangle( const Point & p1,
			    const Color & color1,
			    const Point & p2,
			    const Color & color2,
			    const Point & p3,
			    const Color & color3,
			    unsigned char divisions = 3,
			    unsigned int depth = 0 );

  /** 
   * Draws a triangle with Gouraud-like shaded colors. 
   * 
   * @param x1 
   * @param y1 
   * @param color1 
   * @param x2 
   * @param y2 
   * @param color2 
   * @param x3 
   * @param y3 
   * @param color3 
   * @param divisions 
   * @param depth 
   */
  inline void fillGouraudTriangle( const real x1, const real y1,
				   const Color & color1,
				   const real x2, const real y2, 
				   const Color & color2,
				   const real x3, const real y3,
				   const Color & color3,
				   unsigned char divisions = 3,
				   unsigned int depth = 0 );

  /** 
   * Draws a triangle with a Gouraud-like shaded color according to
   * the current fill color and a brightness value given for each vertex. 
   * @param p1 
   * @param brightness1
   * @param p2 
   * @param brightness2 
   * @param p3 
   * @param brightness3
   * @param divisions number of triangle subdivisions.
   * @param depth The depth of the triangle.
   */
  void fillGouraudTriangle( const Point & p1,
			    const real brightness1,
			    const Point & p2,
			    const real brightness2,
			    const Point & p3,
			    const real brightness3,
			    unsigned char divisions = 3,
			    unsigned int depth = 0 );

  /** 
   * Draws a triangle with a Gouraud-like shaded color according to
   * the current fill color and a brightness value given for each vertex. 
   * 
   * @param x1 
   * @param y1 
   * @param brightness1
   * @param x2 
   * @param y2 
   * @param brightness2
   * @param x3 
   * @param y3 
   * @param brightness3
   * @param divisions
   * @param depth 
   */
  inline void fillGouraudTriangle( const real x1, const real y1,
				   const real brightness1,
				   const real x2, const real y2, 
				   const real brightness2,
				   const real x3, const real y3,
				   const real brightness3,
				   unsigned char divisions = 3,
				   unsigned int depth = 0 );


  /** 
   * Draws a filled triangle.
   * 
   * @param p1 First vertex.
   * @param p2 Second vertex.
   * @param p3 Third vertex.
   * @param depth Depth of the triangle.
   */
  void fillTriangle( const Point & p1,
		     const Point & p2, 
		     const Point & p3, 
		     unsigned int depth = 0 );
  
  /** 
   * Draws a rectangle.
   * 
   * @param x First coordinate of the upper left corner.
   * @param y Second coordinate of the upper left corner.
   * @param width Width of the rectangle.
   * @param height Height of the rectangle.
   * @param depth Depth of the rectangle.
   */
  void drawRectangle( real x, real y, 
		      real width, real height,
		      unsigned int depth = 0 );

  /** 
   * Draws a rectangle filled with the current pen color.
   * 
   * @param x First coordinate of the upper left corner.
   * @param y Second coordinate of the upper left corner.
   * @param width Width of the rectangle.
   * @param height Height of the rectangle.
   * @param depth Depth of the rectangle.
   */
  void fillRectangle( real x, real y,
		      real width, real height,
		      unsigned int depth = 0 );

  /** 
   * Draws a circle.
   * 
   * @param x First coordinate of the circle's center.
   * @param y Second coordinate of the circle's center.
   * @param radius Radius of the circle.
   * @param depth Depth of the circle.
   */
  void drawCircle( real x, real y, real radius,
		   unsigned int depth = 0 );
 
  /** 
   * Draws a circle filled with the current pen color.
   * 
   * @param x First coordinate of the circle's center.
   * @param y Second coordinate of the circle's center.
   * @param radius Radius of the circle.
   * @param depth Depth of the circle.
   */
  void fillCircle( real x, real y, real radius,
		   unsigned int depth = 0);

  /** 
   * Draws an ellipse.
   * 
   * @param x First coordinate of the circle's center.
   * @param y Second coordinate of the circle's center.
   * @param radius Radius of the circle.
   * @param depth Depth of the circle.
   */
  void drawEllipse( real x, real y, 
		    real xRadius, real yRadius,
		    unsigned int depth = 0);
 
  /** 
   * Draws an ellipse filled with the current pen color.
   * 
   * @param x First coordinate of the circle's center.
   * @param y Second coordinate of the circle's center.
   * @param xRadius X axis radius of the ellipse.
   * @param yRadius Y axis radius of the ellipse.
   * @param depth Depth of the circle.
   */
  void fillEllipse( real x, real y, 
		    real xRadius, real yRadius,
		    unsigned int depth = 0);

  /** 
   * Draws a polygonal line.
   * 
   * @param points A vector of points.
   * @param depth The depth of the polyline.
   */
  void drawPolyline( const std::vector<Point> & points,
		     unsigned int depth = 0 );
  
  /** 
   * Draws a closed polygonal line.
   * 
   * @param points A vector of points.
   * @param depth The depth of the polyline.
   */
  void drawClosedPolyline( const std::vector<Point> & points,
			   unsigned int depth = 0 );

  /** 
   * Draws a filled polygon.
   * 
   * @param points A vector of points.
   * @param depth The depth of the polygon.
   */
  void fillPolyline( const std::vector<Point> & points,
		     unsigned int depth = 0 );
    
  /** 
   * Draws a string of text.
   * 
   * @param x The first coordinates of the lower left corner.
   * @param y The second coordinates of the lower left corner.
   * @param text The text. 
   * @param depth The depth of the text.
   */
  void drawText( real x, real y, const char * text, 
		 unsigned int depth = 0 );

  /** 
   * Changes the current font and font size.
   *
   * @param font The name of the font.
   * @param fontSize The new font size.
   * @return The board itself.
   */
  Board & setFont( std::string font, real fontSize );
  
  /** 
   * Changes the font size.
   * 
   * @param fontSize The new font size.
   * @return The board itself.
   */
  Board & setFontSize( real fontSize ); 

  /** 
   * Changes the current pen color.
   * 
   * @param red Red component.
   * @param green Green component.
   * @param blue Blue component.
   * @return The board itself.
   */
  Board & setPenColorRGBi( unsigned char red,
			   unsigned char green,
			   unsigned char blue,
			   unsigned char alpha = 255 );

  /** 
   * Changes the current pen color.
   * 
   * @param red Red
   * @param green 
   * @param blue 
   * @param alpha 
   * @return The board itself.
   */  
  Board & setPenColorRGBf(  real red,
			    real green,
			    real blue, 
			    real alpha = 1.0f );

  /** 
   * Changes the current pen color.
   *
   * In order to use no pen, one may set the pen color to Color::none. 
   * @param color The pen color.
   * @return The board itself.
   */
  Board & setPenColor( const Color & color );
  

  /** 
   * Changes the current fill color.
   * 
   * @param red Red component.
   * @param green Green component.
   * @param blue Blue component.
   * @param alpha The opacity. 
   * @return The board itself.
   */
  Board & setFillColorRGBi( unsigned char red,
			    unsigned char green,
			    unsigned char blue,
			    unsigned char alpha = 255 );
  
  /** 
   * Changes the current fill color.
   * 
   * @param red Red component.
   * @param green Green component.
   * @param blue Blue component.
   * @param alpha The opacity.
   * @return The board itself.
   */
  Board & setFillColorRGBf( real red, real green, real blue, real alpha = 1.0f );

  /** 
   * Changes the current fill color.
   * 
   * In order to use no fill color, one may set this color to Color::none. 
   * @param color The fill color.
   * @return The board itself.
   */
  Board & setFillColor( const Color & color );
  
  /** 
   * Changes the current line thickness (1/72 inche unit).
   * 
   * @param width The new line thickness.
   * @return The board itself.
   */
  Board & setLineWidth( real width );
  
  /** 
   * Set the line cap style. 
   * 
   * @param cap The cap-style which can be Form::ButtCap, 
   * Form::RoundCap or Form::SquareCap.
   * 
   * @return The board itself.
   */  
  inline Board & setLineCap( Form::LineCap cap ); 
 
  /** 
   * Set the line joine style. 
   * 
   * @param cap The join-style which can be Form::MiterJoin, 
   * Form::RoundJoin or Form::BevelJoin.
   * 
   * @return The board itself.
   */  
  inline Board & setLineJoin( Form::LineJoin join );

  /** 
   * Changes the background color of the whole drawing.
   * 
   * @param color A color (may be Color::none).
   */
  void backgroundColor( const Color & color );

  /** 
   * Draws the current drawing's bounding box as a rectangle.
   * 
   * @param depth The depth of the rectangle.
   */
  void drawBoundingBox( unsigned int depth = 0 );

  /** 
   * Save the drawing in an EPS, XFIG of SVG file depending 
   * on the filename extension.
   * 
   * @param filename Path of the file to be created.
   */
  void save( const char * filename ) const; 

  /** 
   * Saves the drawing in an EPS file.
   * 
   * @param filename The EPS file name.
   * @param scale A scale factor to be applied to the while figure before saving.
   */
  void saveEPS( const char * filename, real scale = 1.0f ) const ;

  /** 
   * Saves the drawing in an XFig file.
   * 
   * @param filename The name of the FIG file.
   * @param scale A scale factor to be applied to the while figure before saving.
   */
  void saveFIG( const char * filename, real scale = 1.0f ) const;

  /** 
   * Save the drawing in an SVG file.
   * 
   * @param filename The name of the file.
   * @param scale A scale factor to be applied to the while figure before saving.
   */
  void saveSVG( const char * filename, real scale = 1.0f ) const;

  
 protected:

  /** 
   * Computes the drawing's current bounding box.
   *
   * @return The current bounding box of the drawing.
   */
  Rect computeBoundingBox() const;
  
  Color _penColor;		/**< The current pen color. */
  Color _fillColor;		/**< The current fill color. */
  real _lineWidth;		/**< The current line thickness. */
  Form::LineCap _lineCap;
  Form::LineJoin _lineJoin;
  std::string _font;		/**< The current font. */
  real _fontSize;		/**< The current font size. */
  unsigned int _depth;		/**< The current drawing depth. */
  std::vector< Form* > _shapes; /**< Vector of the shapes. */

  Color _backgroundColor;	/**< The color of the background. */
};

extern const char * XFigPostscriptFontnames[];

inline void
Board::clear( unsigned char red, unsigned char green, unsigned char blue )
{
  clear( Color( red, green, blue ) );
}

inline Board &
Board::setLineCap( Form::LineCap cap )
{
  _lineCap = cap;
  return *this;
}
  
inline Board &
Board::setLineJoin( Form::LineJoin join )
{
  _lineJoin = join;
  return *this;
}

inline void
Board::fillGouraudTriangle( const real x1, const real y1,
			    const Color & color1,
			    const real x2, const real y2, 
			    const Color & color2,
			    const real x3, const real y3,
			    const Color & color3,
			    unsigned char divisions,
			    unsigned int depth )
{
  fillGouraudTriangle( Point( x1, y1 ), color1,
		       Point( x2, y2 ), color2,
		       Point( x3, y3 ), color3,
		       divisions, depth );		       
}

void
Board::fillGouraudTriangle( const real x1, const real y1,
			    const real brightness1,
			    const real x2, const real y2, 
			    const real brightness2,
			    const real x3, const real y3,
			    const real brightness3,
			    unsigned char divisions,
			    unsigned int depth )
{
  fillGouraudTriangle( Point( x1, y1 ), brightness1,
		       Point( x2, y2 ), brightness2,
		       Point( x3, y3 ), brightness3,
		       divisions, depth );
}

#endif
