#include "board.h"
#include "point.h"
#include "rect.h"
#include "shapes.h"
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <map>
#include <algorithm>
#include <cstdio>
#include <algorithm>

const char * XFigPostscriptFontnames[] = {
  "Times-Roman",
  "Times-Italic",
  "Times-Bold",
  "Times-Bold-Italic",
  "AvantGarde-Book",
  "AvantGarde-Book-Oblique",
  "AvantGarde-Demi",
  "AvantGarde-Demi-Oblique",
  "Bookman-Light",
  "Bookman-Light-Italic",
  "Bookman-Demi",
  "Bookman-Demi-Italic",
  "Courier",
  "Courier-Oblique",
  "Courier-Bold",
  "Courier-Bold-Oblique",
  "Helvetica",
  "Helvetica-Oblique",
  "Helvetica-Bold",
  "Helvetica-Bold-Oblique",
  "Helvetica-Narrow",
  "Helvetica-Narrow-Oblique",
  "Helvetica-Narrow-Bold",
  "Helvetica-Narrow-Bold-Oblique",
  "New-Century-Schoolbook-Roman",
  "New-Century-Schoolbook-Italic",
  "New-Century-Schoolbook-Bold",
  "New-Century-Schoolbook-Bold-Italic",
  "Palatino-Roman",
  "Palatino-Italic",
  "Palatino-Bold",
  "Palatino-Bold-Italic",
  "Symbol",
  "Zapf-Chancery-Medium-Italic",
  "Zapf-Dingbats",
  NULL
};

Board::Board( const Color & backgroundColor )
  :_penColor( Color::black ),
   _fillColor( Color::none ),
   _lineWidth(0.5),
   _lineCap( Form::ButtCap ),
   _lineJoin( Form::MiterJoin ),
   _font( std::string( "Times-Roman" ) ),
   _fontSize( 11.0 ),
   _depth( std::numeric_limits<unsigned int>::max() - 1 ),
   _backgroundColor( backgroundColor )
{ 
}

Board::~Board() {
  std::vector< Form* >::iterator i = _shapes.begin();
  std::vector< Form* >::iterator end = _shapes.end();
  while ( i != end ) {
    delete *i;
    *i = 0;
    ++i;
  }
}

void
Board::clear( const Color & color )
{
  _shapes.clear();
  _depth = std::numeric_limits<unsigned int>::max() - 1;
  _penColor.setRGBi( 0, 0, 0);
  _fillColor = Color::none;
  _backgroundColor = color;
}

Board &
Board::setPenColorRGBi( unsigned char red,
			unsigned char green,
			unsigned char blue, 
			unsigned char alpha )
{
  _penColor.setRGBi( red, green, blue, alpha );
  return *this;
}

Board &
Board::setPenColorRGBf(  real red,
			 real green,
			 real blue, 
			 real alpha )
{
  if ( red > 1.0f ) red = 1.0f;
  if ( red < 0.0f ) red = 0.0f;
  if ( green > 1.0f ) green = 1.0f;
  if ( green < 0.0f ) green = 0.0f;
  if ( blue > 1.0f ) blue = 1.0f;
  if ( blue > 1.0f ) red = 1.0;
  return *this;
}

Board &
Board::setPenColor( const Color & color )
{
  _penColor = color;
  return *this;
}

Board &
Board::setFillColorRGBi( unsigned char red,
			 unsigned char green,
			 unsigned char blue,
			 unsigned char alpha )
{
  _fillColor.setRGBi( red, green, blue, alpha );
  return *this;
}

Board &
Board::setFillColorRGBf( real red, real green, real blue, real alpha )
{
  _fillColor.setRGBf( red, green, blue, alpha );
  return *this;
}

Board &
Board::setFillColor( const Color & color )
{
  _fillColor = color;
  return *this;
}

Board &
Board::setLineWidth( real width )
{
  _lineWidth = width;
  return *this;
}

Board &
Board::setFont( std::string font, real fontSize )
{
  _font = font;
  _fontSize = fontSize;
  return *this;
}

Board &
Board::setFontSize( real fontSize )
{
  _fontSize = fontSize;
  return *this;
}

void
Board::backgroundColor( const Color & color )
{
  _backgroundColor = color;
}

void
Board::drawLine( real x1, real y1, real x2, real y2, 
		   unsigned int depth /* = 0 */  )
{
  if ( depth ) 
    _shapes.push_back( new Line( x1, y1, x2, y2, _penColor, _lineWidth, _lineCap, _lineJoin, depth ) );
  else
    _shapes.push_back( new Line( x1, y1, x2, y2, _penColor, _lineWidth, _lineCap, _lineJoin, _depth-- ) );
}

void
Board::drawArrow( real x1, real y1, real x2, real y2, 
		  bool filled /* = false */,
		  unsigned int depth /* = 0 */  )
{
  if ( depth )
    _shapes.push_back( new Arrow( x1, y1, x2, y2,
				  _penColor, filled ? _penColor : Color::none,
				  _lineWidth, _lineCap, _lineJoin, depth ) );
  else
    _shapes.push_back( new Arrow( x1, y1, x2, y2,
				  _penColor, filled ? _penColor : Color::none,
				  _lineWidth, _lineCap, _lineJoin, _depth-- ) );
}

void
Board::drawRectangle( real x, real y, 
		      real width, real height,
		      unsigned int depth /* = 0 */ )
{

  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Rectangle( x, y, width, height, 
				    _penColor, _fillColor,
				    _lineWidth, _lineCap, _lineJoin, d ) );
}

void
Board::fillRectangle( real x, real y,
		      real width, real height,
		      unsigned int depth /* = 0 */ )
{

  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Rectangle( x, y, width, height, 
				    Color::none, _penColor,
				    0.0f, _lineCap, _lineJoin,
				    d ) );

}

void
Board::drawCircle( real x, real y, real radius,
		   unsigned int depth /* = 0 */  )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Circle( x, y, radius, 
				 _penColor, _fillColor,
				 _lineWidth, _lineCap, _lineJoin,
				 d ) );
}

void
Board::fillCircle( real x, real y, real radius,
		   unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Circle( x, y, radius, 
				 Color::none, _penColor,
				 0.0f, _lineCap, _lineJoin,
				 d ) );
}

void
Board::drawEllipse( real x, real y,
		    real xRadius, real yRadius,
		    unsigned int depth /* = 0 */  )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Ellipse( x, y, xRadius, yRadius,
				  _penColor, _fillColor,
				  _lineWidth, _lineCap, _lineJoin,
				  d ) );
}

void
Board::fillEllipse( real x, real y, 
		    real xRadius, real yRadius,
		    unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Ellipse( x, y, xRadius, yRadius,
				  Color::none, _penColor,
				  0.0f, _lineCap, _lineJoin,
				  d ) );
}

void
Board::drawPolyline( const std::vector<Point> & points,
		       unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Polyline( points, false, _penColor, _fillColor,
				   _lineWidth, _lineCap, _lineJoin,
				   d ) );
}

void
Board::drawClosedPolyline( const std::vector<Point> & points,
			     unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Polyline( points, true, _penColor, _fillColor,
				   _lineWidth, _lineCap, _lineJoin,
				   d ) );
}

void
Board::fillPolyline( const std::vector<Point> & points,
		       unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Polyline( points, true, Color::none, _penColor,
				   0.0f, _lineCap, _lineJoin,
				   d ) );
}

void
Board::drawTriangle( real x1, real y1, 
		     real x2, real y2, 
		     real x3, real y3, 
		     unsigned int depth )
{
  unsigned int d = depth ? depth : _depth--;
  std::vector<Point> points;
  points.push_back( Point( x1, y1 ) );
  points.push_back( Point( x2, y2 ) );
  points.push_back( Point( x3, y3 ) );
  _shapes.push_back( new Polyline( points, true, _penColor, _fillColor,
				   _lineWidth, _lineCap, _lineJoin,
				   d ) );
}

void
Board::drawTriangle( const Point & p1,
		     const Point & p2, 
		     const Point & p3, 
		     unsigned int depth )
{
  unsigned int d = depth ? depth : _depth--;
  std::vector<Point> points;
  points.push_back( p1 );
  points.push_back( p2 );
  points.push_back( p3 );
  _shapes.push_back( new Polyline( points, true, _penColor, _fillColor,
				   _lineWidth, _lineCap, _lineJoin,
				   d ) );  
}

void
Board::fillTriangle( real x1, real y1, 
		     real x2, real y2, 
		     real x3, real y3, 
		     unsigned int depth )
{
  unsigned int d = depth ? depth : _depth--;
  std::vector<Point> points;
  points.push_back( Point( x1, y1 ) );
  points.push_back( Point( x2, y2 ) );
  points.push_back( Point( x3, y3 ) );
  _shapes.push_back( new Polyline( points, true, Color::none, _penColor,
				   0.0f, _lineCap, _lineJoin,
				   d ) );
}

void
Board::fillTriangle( const Point & p1,
		     const Point & p2, 
		     const Point & p3, 
		     unsigned int depth )
{
  unsigned int d = depth ? depth : _depth--;
  std::vector<Point> points;
  points.push_back( p1 );
  points.push_back( p2 );
  points.push_back( p3 );
  _shapes.push_back( new Polyline( points, true, Color::none, _penColor,
				   0.0f, _lineCap, _lineJoin,
				   d ) );  
}

void
Board::fillGouraudTriangle( const Point & p1,
			    const Color & color1,
			    const Point & p2,
			    const Color & color2,
			    const Point & p3,
			    const Color & color3,
			    unsigned char divisions,
			    unsigned int depth )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new GouraudTriangle( p1, color1, p2, color2, p3, color3,
					  divisions, d ) );
}

void
Board::fillGouraudTriangle( const Point & p1,
			    const real brightness1,
			    const Point & p2,
			    const real brightness2,
			    const Point & p3,
			    const real brightness3,
			    unsigned char divisions,
			    unsigned int depth )
{
  Color color1( _penColor );
  Color color2( _penColor );
  Color color3( _penColor );
  color1.red( static_cast<unsigned char>( std::min( 255.0, color1.red() * brightness1 ) ) );
  color1.green( static_cast<unsigned char>( std::min( 255.0, color1.green() * brightness1 ) ) );
  color1.blue( static_cast<unsigned char>( std::min( 255.0, color1.blue() * brightness1 ) ) );
  color2.red( static_cast<unsigned char>( std::min( 255.0, color2.red() * brightness2 ) ) );
  color2.green( static_cast<unsigned char>( std::min( 255.0, color2.green() * brightness2 ) ) );
  color2.blue( static_cast<unsigned char>( std::min( 255.0, color2.blue() * brightness2 ) ) );
  color3.red( static_cast<unsigned char>( std::min( 255.0, color3.red() * brightness3 ) ) );
  color3.green( static_cast<unsigned char>( std::min( 255.0, color3.green() * brightness3 ) ) );
  color3.blue( static_cast<unsigned char>( std::min( 255.0, color3.blue() * brightness3 ) ) );
  fillGouraudTriangle( p1, color1, p2, color2, p3, color3, divisions, depth );
}

void
Board::drawText( real x, real y, const char * text,
		 unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  _shapes.push_back( new Text( x, y, text, _font, _fontSize, _penColor, d ) );
}

void
Board::drawBoundingBox( unsigned int depth /* = 0 */ )
{
  unsigned int d = depth ? depth : _depth--;
  Rect bbox = computeBoundingBox();
  _shapes.push_back( new Rectangle( bbox.left,
				    bbox.top,
				    bbox.width,
				    bbox.height,
				    _penColor, _fillColor,
				    _lineWidth, _lineCap, _lineJoin,
				    d ) );
}

Rect
Board::computeBoundingBox() const
{
  Rect result;
  std::vector< Form* >::const_iterator i = _shapes.begin();
  std::vector< Form* >::const_iterator end = _shapes.end();
  if ( i == end ) return result;
  result = (*i)->boundingBox();
  ++i;
  while ( i != end ) { 
    result = result || (*i)->boundingBox();
    ++i;
  }
  return result;
}
  
void
Board::saveEPS( const char * filename, real scale ) const
{
  std::ofstream file( filename );
  
  Rect bbox = computeBoundingBox();
  TransformEPS transform;
  transform.setBoundingBox( bbox );
  if ( scale != 1.0 ) transform.furtherScale( scale );

  file << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
  file << "%%Title: " << filename << std::endl;
  file << "%%Creator: Board library (C)2007 Sebastien Fourey" << std::endl;
  time_t t = time(0);
  file << "%%CreationDate: " << ctime(&t);
  file << "%%BoundingBox: " << std::setprecision( 8 )
       << transform.mapX( bbox.left ) << " "
       << transform.mapY( bbox.top - bbox.height ) << " "
       << transform.mapX( bbox.left + bbox.width ) << " "
       << transform.mapY( bbox.top ) << std::endl;
  file << "%Magnification: 1.0000" << std::endl;
  file << "%%EndComments" << std::endl;

  file << std::endl;
  file << "/cp {closepath} bind def" << std::endl;
  file << "/ef {eofill} bind def" << std::endl;
  file << "/gr {grestore} bind def" << std::endl;
  file << "/gs {gsave} bind def" << std::endl;
  file << "/sa {save} bind def" << std::endl;
  file << "/rs {restore} bind def" << std::endl;
  file << "/l {lineto} bind def" << std::endl;
  file << "/m {moveto} bind def" << std::endl;
  file << "/rm {rmoveto} bind def" << std::endl;
  file << "/n {newpath} bind def" << std::endl;
  file << "/s {stroke} bind def" << std::endl;
  file << "/sh {show} bind def" << std::endl;
  file << "/slc {setlinecap} bind def" << std::endl;
  file << "/slj {setlinejoin} bind def" << std::endl;
  file << "/slw {setlinewidth} bind def" << std::endl;
  file << "/srgb {setrgbcolor} bind def" << std::endl;
  file << "/rot {rotate} bind def" << std::endl;
  file << "/sc {scale} bind def" << std::endl;
  file << "/sd {setdash} bind def" << std::endl;
  file << "/ff {findfont} bind def" << std::endl;
  file << "/sf {setfont} bind def" << std::endl;
  file << "/scf {scalefont} bind def" << std::endl;
  file << "/sw {stringwidth} bind def" << std::endl;
  file << "/tr {translate} bind def" << std::endl;
  file << " 0.5 setlinewidth" << std::endl;

  file << "newpath " 
       << transform.mapX( bbox.left ) 
       << " " << transform.mapY( bbox.top ) << " m " 
       << transform.mapX( bbox.left )
       << " " << transform.mapY( bbox.top - bbox.height ) << " l " 
       << transform.mapX( bbox.left + bbox.width  ) << " "
       << transform.mapY( bbox.top - bbox.height ) << " l " 
       << transform.mapX( bbox.left + bbox.width ) << " "
       << transform.mapY( bbox.top ) << " l " 
       << " cp clip " << std::endl;

  
  // Draw the background color if needed.
  if ( _backgroundColor != Color::none ) { 
    Rectangle r( bbox, Color::none, _backgroundColor, 0.0f );
    r.flushPostscript( file, transform );
  }

  // Draw the shapes
  std::vector< Form* > shapes = _shapes;

  stable_sort( shapes.begin(), shapes.end(), shapeGreaterDepth );
  std::vector< Form* >::const_iterator i = shapes.begin();
  std::vector< Form* >::const_iterator end = shapes.end();

  Color color( 0, 0, 0 );
  real lineWidth = 1.0f;
  Form::LineCap lineCap = Form::ButtCap;
  Form::LineJoin lineJoin = Form::MiterJoin;

  while ( i != end ) {
    if ( (*i)->penColor != color ) { 
      color = (*i)->penColor;
      color.flushPostscript( file );
    }
    if ( (*i)->lineWidth != lineWidth ) {
      lineWidth = (*i)->lineWidth;
      file << lineWidth << " slw ";
    }
    if ( (*i)->lineCap != lineCap ) {
      lineCap = (*i)->lineCap;
      file << lineCap << " slc ";
    }
    if ( (*i)->lineJoin != lineJoin ) {
      lineJoin = (*i)->lineJoin;
      file << lineJoin << " slj ";
    }
    (*i)->flushPostscript( file, transform );
    ++i;
  }
  file << "showpage" << std::endl;
  file << "%%Trailer" << std::endl;
  file << "%EOF" << std::endl;
  file.close();
}

void
Board::saveFIG( const char * filename, real scale ) const
{
  std::ofstream file( filename );
  TransformFIG transform;
  Rect bbox = computeBoundingBox();
  transform.setBoundingBox( bbox );
  transform.setDepthRange( _shapes );
  if ( scale != 1.0 ) transform.furtherScale( scale );

  file << "#FIG 3.2  Produced by the Board library (C)2007 Sebastien Fourey\n";
  file << "Portrait\n";
  file << "Center\n";
  file << "Metric\n";
  file << "A4\n";
  file << "100.00\n";
  file << "Single\n";
  file << "-2\n";
  file << "1200 2\n";

  std::map<Color,int> colormap;
  int maxColor = 32;


  colormap[Color(0,0,0)] = 0; 
  colormap[Color(0,0,255)] = 1; 
  colormap[Color(0,255,0)] = 2; 
  colormap[Color(0,255,255)] = 0; 
  colormap[Color(255,0,0)] = 4; 
  colormap[Color(255,0,255)] = 0; 
  colormap[Color(255,255,0)] = 6; 
  colormap[Color(255,255,255)] = 7;


  std::vector< Form* > shapes = _shapes;
  stable_sort( shapes.begin(), shapes.end(), shapeGreaterDepth );
  std::vector< Form* >::const_iterator i = shapes.begin();
  std::vector< Form* >::const_iterator end = shapes.end();
  while ( i != end ) { 
    if ( colormap.find( (*i)->penColor ) == colormap.end() 
	 && (*i)->penColor.valid() )
      colormap[ (*i)->penColor ] = maxColor++;
    if ( colormap.find( (*i)->fillColor ) == colormap.end()
	 && (*i)->fillColor.valid() )
      colormap[ (*i)->fillColor ] = maxColor++;
    ++i;
  }

  if ( colormap.find( _backgroundColor ) == colormap.end()
       && _backgroundColor.valid() )
    colormap[ _backgroundColor ] = maxColor++;
  

  // Write the colormap
  std::map<Color,int>::const_iterator iColormap = colormap.begin();
  std::map<Color,int>::const_iterator endColormap = colormap.end();
  char colorString[255];
  while ( iColormap != endColormap ) {
    sprintf( colorString,"0 %d #%02x%02x%02x\n",
	     iColormap->second,
	     iColormap->first.red(),
	     iColormap->first.green(),
	     iColormap->first.blue() );
    if ( iColormap->second >= 32 ) file << colorString;
    ++iColormap;
  }

  // Draw the background color if needed.
  if ( _backgroundColor != Color::none ) { 
    Rectangle r( bbox, Color::none, _backgroundColor, 0.0f );
    r.flushFIG( file, transform, colormap );
  }
  // Draw the shapes.
  i = shapes.begin();
  while ( i != end ) {
    (*i)->flushFIG( file, transform, colormap );
    ++i;
  }  
  file.close();
}

void
Board::saveSVG( const char * filename, real scale ) const
{
  real ppmm = 720/254.0f;
  std::ofstream file( filename );
  TransformSVG transform;
  Rect bbox = computeBoundingBox();
  transform.setBoundingBox( bbox );
  if ( scale != 1.0 ) transform.furtherScale( scale );

  file << "<?xml version=\"1.0\" standalone=\"no\"?>" << std::endl;
  file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
  file << " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;
  file << "<svg width=\"210mm\" height=\"297mm\" " << std::endl;
  file << "     viewBox=\"0 0 " << 210*ppmm << " " << 297*ppmm << "\" " << std::endl;
  file << "     xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << std::endl;
  file << "<desc>" << filename << ", created with the Board library (C) 2007 Sebastien Fourey" << "</desc>" << std::endl;
  
  // Draw the background color if needed.
  if ( _backgroundColor != Color::none ) { 
    Rectangle r( bbox, Color::none, _backgroundColor, 0.0 );
    r.flushSVG( file, transform );
  }
  
  // Draw the shapes.
  std::vector< Form* > shapes = _shapes;
  stable_sort( shapes.begin(), shapes.end(), shapeGreaterDepth );
  std::vector< Form* >::const_iterator i = shapes.begin();
  std::vector< Form* >::const_iterator end = shapes.end();
  while ( i != end ) {
    (*i)->flushSVG( file, transform );
    ++i;
  }  
  file << "</svg>" << std::endl;
  file.close();
}

void
Board::save( const char * filename ) const 
{
  const char * extension = filename + strlen( filename );
  while ( extension > filename && *extension != '.' ) 
    --extension;
  if ( !(strcmp( extension, ".eps" )) || !(strcmp( extension, ".EPS" )) ) {
    saveEPS( filename );
    return;
  }
  if ( !(strcmp( extension, ".fig" )) || !(strcmp( extension, ".FIG" )) ) {
    saveFIG( filename );
    return;
  }
  if ( !(strcmp( extension, ".svg" )) || !(strcmp( extension, ".SVG" )) ) {
    saveSVG( filename );
    return;
  }
}

/**
 * @example examples/arithmetic.cc
 * @example examples/example1.cc
 * @example examples/example2.cc
 */

/**
 * @mainpage Board - A C++ library for simple Postscript, SVG, and XFig drawings.
 *
 * <img align=left src="http://www.greyc.ensicaen.fr/~seb/images/board.png"> (C) 2007 Sébastien Fourey - GREYC ENSICAEN 
 *
 * @section intro_sec Introduction
 *
 * The board library allows simple drawings in:
 * <ul>
 *  <li>Encapsulated Postcript files (EPS) ;
 *  <li>XFig files (FIG) ;
 *  <li>Scalable Vector Graphics files (SVG).
 * </ul>
 *
 * The main class of the library is the #BoardLib#Board class. It is intented to be as simple as possible
 * so that it can be used quickly in programs to generate the kind of figure one would rather
 * not draw by hand, but which can be easily drawn by a computer (C++) program.
 *
 * @section examples_sec Code examples
 *
 * <ul>
 * <li>Se the "Examples" tab above.
 * </ul>
 *
 * @section links_sec Links
 *
 * <ul>
 * <li>Visit the <a href="http://www.greyc.ensicaen.fr/~seb/board/">Board homepage</a>.</li>
 * </ul>
 */
