#ifndef _BOARD_RECT_H_
#define _BOARD_RECT_H_

#include "../Global/global_precision.h"

/**
 * The Rectangle structure.
 * @brief Struct representing a rectangle on the plane. 
 */
struct Rect {
  real left;			/**< Coordinate of the left side. */
  real top;			/**< Coordinate of the upper side. */
  real width;			/**< Width of the rectangle. */
  real height;		/**< Height of the rectangle. */
  
  /** 
   * Rect constructor.
   * 
   * @param left 
   * @param top 
   * @param width 
   * @param height 
   * 
   * @return 
   */
  Rect( real left = 0.0f, real top = 0.0f, real width = 0.0f, real height = 0.0f )
    :left( left ), top( top ), width( width ), height( height ) { } 
};

/** 
 * Computes the bounding box of two bounding boxes.
 * 
 * @param rectA A first rectangle.
 * @param rectB A second rectangle.
 * 
 * @return The smallest rectangle that contains both rectA and rectB.
 */
Rect operator||( const Rect & rectA, const Rect & rectB );

#endif // _RECT_H_

/*-----------------------------------------------------------------------------+
 '$Id: rect.h,v 1.3 2008/11/17 19:23:22 niceno Exp $'/
+-----------------------------------------------------------------------------*/
