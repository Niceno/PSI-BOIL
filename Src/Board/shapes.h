#ifndef _BOARD_SHAPES_H_
#define _BOARD_SHAPES_H_

#include "point.h"
#include "rect.h"
#include "color.h"
#include "transforms.h"
#include <string>
#include <vector>
#include <iostream>
#include <map>

/**
 * Form structure.
 * @brief Abstract structure for a 2D shape.
 */
struct Form {

  enum LineCap { ButtCap = 0, RoundCap, SquareCap };
  enum LineJoin { MiterJoin = 0, RoundJoin, BevelJoin };

  unsigned int depth;		/**< The depth of the shape. */
  Color penColor;		/**< The color of the shape. */
  Color fillColor;		/**< The color of the shape. */
  real lineWidth;		/**< The line thickness. */
  LineCap lineCap;		/**< The linecap attribute. (The way line terminates.) */
  LineJoin lineJoin;		/**< The linejoin attribute. (The shape of line junctions.) */

  /** 
   * Form constructor.
   * 
   * @param penColor The pen color of the shape.
   * @param fillColor The fill color of the shape.
   * @param lineWidth The line thickness.
   * @param depth The depth of the shape.
   */
  Form( Color penColor, Color fillColor,
	 real lineWidth, const LineCap cap, const LineJoin join,
	 unsigned int depth )
    : depth( depth ), penColor( penColor ), fillColor( fillColor ), 
      lineWidth( lineWidth ), lineCap( cap ), lineJoin( join ) { }

  /** 
   * Form destructor.
   */
  virtual ~Form() { }
  
  /** 
   * Checks whether a shape is filled with a color or not.
   * 
   * @return true if the shape is filled.
   */
  inline bool filled() const { return fillColor != Color::none; }

  /** 
   * Writes the EPS code of the shape in a stream according
   * to a transform.
   * 
   * @param stream The output stream.
   * @param transform A 2D transform to be applied.
   */
  virtual void flushPostscript( std::ostream & stream,
				const TransformEPS & transform ) const = 0;

  /** 
   * Writes the FIG code of the shape in a stream according
   * to a transform.
   * 
   * @param stream The output stream.
   * @param transform A 2D transform to be applied.
   */
  virtual void flushFIG( std::ostream & stream,
			 const TransformFIG & transform,
			 std::map<Color,int> & colormap ) const = 0;

  /** 
   * Writes the SVG code of the shape in a stream according
   * to a transform.
   * 
   * @param stream The output stream.
   * @param transform A 2D transform to be applied.
   */
  virtual void flushSVG( std::ostream & stream,
			 const TransformSVG & transform ) const = 0;
  
  /** 
   * Returns the bounding box of the figure.
   *
   * @return The rectangle of the bounding box.
   */
  virtual Rect boundingBox() const = 0;

protected:

  /** 
   * Return a string of the svg properties lineWidth, opacity, penColor, fillColor,
   * lineCap, and lineJoin.
   * 
   * @return A string of the properties suitable for inclusion in an svg tag.
   */
  std::string svgProperties( const TransformSVG & transform ) const;
};

/**
 * The line structure.
 * @brief A line between two points.
 */
struct Line : public Form { 
  real x1;			/**< First coordinate of the start point. */
  real y1;			/**< Second coordinate of the start point. */
  real x2; 			/**< First coordinate of the end point. */
  real y2;			/**< Second coordinate of the end point. */
  
  /** 
   * Constructs a line.
   * 
   * @param x1 First coordinate of the start point.
   * @param y1 Second coordinate of the start point.
   * @param x2 First coordinate of the end point.
   * @param y2 Second coordinate of the end point.
   * @param color The color of the line.
   * @param lineWidth The line thickness.
   * @param depth The depth of the line.
   */
  Line( real x1, real y1, real x2, real y2, 
	Color color, 
	real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	unsigned int depth = 0 )
    : Form( color, Color::none, lineWidth, cap, join, depth ),
      x1( x1 ), y1( y1 ), x2( x2 ), y2( y2 ) { }

  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;
  
  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};

/**
 * The arrow structure.
 * @brief A line between two points with an arrow at one extremity.
 */
struct Arrow : public Line { 

  /** 
   * Constructs an arrow.
   * 
   * @param x1 First coordinate of the start point.
   * @param y1 Second coordinate of the start point.
   * @param x2 First coordinate of the end point.
   * @param y2 Second coordinate of the end point.
   * @param penColor The color of the line.
   * @param fillColor The fill color of the sharp end.
   * @param lineWidth The line thickness.
   * @param depth The depth of the line.
   */
  Arrow( real x1, real y1, real x2, real y2,
	 Color penColor, Color fillColor,
	 real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	 unsigned int depth = 0 )
    : Line( x1, y1, x2, y2, penColor, lineWidth, cap, join, depth ) {
    Form::fillColor = fillColor;
  }
  
  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;
  
  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;
  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

};

/**
 * The polyline structure.
 * @brief A polygonal line described by a series of 2D points.
 */
struct Polyline : public Form { 
  std::vector<Point> points;
  bool closed;

  Polyline( const std::vector<Point> & points, 
	    bool closed,
	    Color penColor, Color fillColor,
	    real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	    unsigned int depth = 0 )
    : Form( penColor, fillColor, lineWidth, cap, join, depth ),
      points( points ), closed( closed ) { }
	    
  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};


/**
 * The GouraudTriangle structure.
 * @brief A triangle with shaded filling according to colors given for each vertex. 
 */
struct GouraudTriangle : public Polyline {
  Color color0;
  Color color1;
  Color color2;
  int subdivisions;

  GouraudTriangle( const Point & p0, const Color & color0,
		   const Point & p1, const Color & color1,
		   const Point & p2, const Color & color2,
		   int subdivisions,
		   unsigned int depth = 0 );

  GouraudTriangle( const Point & p0, real brightness0,
		   const Point & p1, real brightness1,
		   const Point & p2, real brightness2,
		   const Color & fillColor,
		   int subdivisions,
		   unsigned int depth = 0 );

  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  /** 
   * Sends the Triangle to a FIG file format stream.
   * <p><b>Warning!</b> Because shading would generally require
   * more colors in the colormap than allowed by the FIG file format, 
   * rendering a Gouraud triangle in an XFig file is the same as rendering
   * a simple triangle whose filling color is the average of the vertex colors.
   * 
   * @param stream 
   * @param transform 
   * @param Color 
   * @param colormap 
   */
  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

};

/**
 * The rectangle structure.
 * @brief A rectangle.
 */
struct Rectangle : public Form {
  real x;
  real y;
  real width;
  real height;

  Rectangle( real x, real y, real width, real height,
	     Color penColor, Color fillColor,
	     real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	     unsigned int depth = 0 )
    : Form( penColor, fillColor, lineWidth, cap, join, depth ),
      x( x ), y( y ), width( width ), height( height ) { }

  Rectangle( const Rect & rect,
	     Color penColor, Color fillColor,
	     real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	     unsigned int depth = 0 )
    : Form( penColor, fillColor, lineWidth, cap, join, depth ),
      x( rect.left ), y( rect.top ), width( rect.width ), height( rect.height ) { }
  
  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};

/**
 * The circle structure.
 * @brief A circle.
 */
struct Circle : public Form {
  real x;
  real y;
  real radius;
  
  Circle( real x, real y, real radius, 
	  Color penColor, Color fillColor,
	  real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	  unsigned int depth = 0 )
    : Form( penColor, fillColor, lineWidth, cap, join, depth ),
      x( x ), y( y ), radius( radius ) { }

  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};

/**
 * The ellipse structure.
 * @brief An ellipse.
 */
struct Ellipse : public Form {
  real x;
  real y;
  real xRadius; 
  real yRadius;
  
  Ellipse( real x, real y, 
	   real xRadius, real yRadius, 
	   Color penColor, Color fillColor,
	   real lineWidth, const LineCap cap = ButtCap, const LineJoin join = MiterJoin,
	   unsigned int depth = 0 )
    : Form( penColor, fillColor, lineWidth, cap, join, depth ),
      x( x ), y( y ), xRadius( xRadius ), yRadius( yRadius ) { }

  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};

/**
 * The text structure.
 * @brief A piece of text.
 */
struct Text : public Form {
  real x;
  real y;
  std::string text;
  std::string font;
  real size;
  
  Text( real x, real y,
	const std::string & text, const std::string & font, real size,
	Color color, unsigned int depth )
    : Form( color, Color::none, 1.0, ButtCap, MiterJoin, depth ),
      x( x ), y( y ), text( text ), font( font ), size( size ) { }
  
  void flushPostscript( std::ostream & stream,
			const TransformEPS & transform ) const;

  void flushFIG( std::ostream & stream,
		 const TransformFIG & transform,
		 std::map<Color,int> & colormap ) const;

  void flushSVG( std::ostream & stream,
		 const TransformSVG & transform ) const;

  Rect boundingBox() const;
};

/** 
 * Compares two shapes according to their depths.
 * 
 * @param s1 A pointer to a first shape.
 * @param s2 A pointer to a second shape.
 * 
 * @return 
 */
bool shapeGreaterDepth( const Form *s1, const Form *s2 );

#endif /* _SHAPE_H_ */

/*-----------------------------------------------------------------------------+
 '$Id: shapes.h,v 1.2 2008/10/21 11:54:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
