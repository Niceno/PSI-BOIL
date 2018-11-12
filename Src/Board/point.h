#ifndef _BOARD_POINT_H_
#define _BOARD_POINT_H_

/**
 * The Point structure.
 * @brief Struct representing a 2D point. 
 */
struct Point {
  real x;			/**< The point's first coordinate */
  real y;			/**< The point's second coordinate */
  /** 
   * Point constructor.
   * 
   * @param x The point's first coordinate.
   * @param y The point's second coordinate.
   */
  Point( real x, real y):x(x),y(y) { } 
};

#endif // _POINT_H_
