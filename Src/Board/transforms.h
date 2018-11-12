#ifndef _BOARD_TRANSFORMS_H_
#define _BOARD_TRANSFORMS_H_

#include <limits>
#include <vector>

struct Rect;
struct Form;

/**
 * The base class for transforms.
 * @brief 
 */
struct Transform {
public:
  Transform():_scale(1.0), _deltaX(0.0), _deltaY(0.0) { }
  virtual ~Transform() { };
  virtual real mapX( real x ) const;
  virtual real mapY( real y ) const = 0;
  virtual void apply( real & x, real & y ) const;
  virtual real scale( real x ) const;
  virtual void furtherScale( real x );
  virtual real rounded( real x ) const;
  virtual void setBoundingBox( const Rect & rect ) = 0;
protected:
  real _scale;
  real _deltaX;
  real _deltaY;
};

/**
 * The TransformEPS structure.
 * @brief Structure representing a scaling and translation
 * suitable for an EPS output.
 */
struct TransformEPS : public Transform {
public:
  real mapY( real y ) const;
  void setBoundingBox( const Rect & rect );
};

/**
 * The TransformFIG structure.
 * @brief Structure representing a scaling and translation
 * suitable for an XFig output.
 */
struct TransformFIG : public Transform {
public:
  TransformFIG():maxDepth(std::numeric_limits<unsigned int>::max()),minDepth(0) { }
  real rounded( real x ) const;
  real mapY( real y ) const;
  int mapWidth( real width ) const; 
  void setBoundingBox( const Rect & rect );
  void setDepthRange( const std::vector<Form*> & shapes );
  unsigned int mapDepth( unsigned int depth ) const;
private:
  unsigned int maxDepth;
  unsigned int minDepth;
};

/**
 * The TransformSVG structure.
 * @brief Structure representing a scaling and translation
 * suitable for an SVG output.
 */
struct TransformSVG : public Transform {
public:
  real rounded( real x ) const;
  real mapY( real y ) const;
  real mapWidth( real width ) const; 
  void setBoundingBox( const Rect & rect );
};

#endif /* _TRANSFORMS_H_ */
