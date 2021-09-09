// https://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle
#include "custom.h"

real setup_section(real h, real r) // returns the positive root of intersection of line y = h with circle centered at the origin and radius r
{
    assert(r >= 0); // assume r is positive, leads to some simplifications in the formula below (can factor out r from the square root)
    return (h < r)? sqrt(r * r - h * h) : 0; // http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+%3D+h
}

real setup_g(real x, real h, real r) // indefinite integral of circle segment
{
    return .5 * (sqrt(1 - x * x / (r * r)) * x * r + r * r * asin(x / r) - 2 * h * x); // http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+-+h
}

real setup_area(real x0, real x1, real h, real r) // area of intersection of an infinitely tall box with left edge at x0, right edge at x1, bottom edge at h and top edge at infinity, with circle centered at the origin with radius r
{
    if(x0 > x1)
        std::swap(x0, x1); // this must be sorted otherwise we get negative area
    real s = setup_section(h, r);
    return setup_g(std::max(-s, std::min(s, x1)), h, r) - setup_g(std::max(-s, std::min(s, x0)), h, r); // integrate the area
}

real setup_area(real x0, real x1, real y0, real y1, real r) // area of the intersection of a finite box with a circle centered at the origin with radius r
{
    if(y0 > y1)
        std::swap(y0, y1); // this will simplify the reasoning
    if(y0 < 0) {
        if(y1 < 0)
            return setup_area(x0, x1, -y0, -y1, r); // the box is completely under, just flip it above and try again
        else
            return setup_area(x0, x1, 0, -y0, r) + setup_area(x0, x1, 0, y1, r); // the box is both above and below, divide it to two boxes and go again
    } else {
        assert(y1 >= 0); // y0 >= 0, which means that y1 >= 0 also (y1 >= y0) because of the swap at the beginning
        return setup_area(x0, x1, y0, r) - setup_area(x0, x1, y1, r); // area of the lower box minus area of the higher box
    }
}

real setup_area(real x0, real x1, real y0, real y1, real cx, real cy, real r) // area of the intersection of a general box with a general circle
{
    x0 -= cx; x1 -= cx;
    y0 -= cy; y1 -= cy;
    // get rid of the circle center

    return setup_area(x0, x1, y0, y1, r);
}

namespace boil {
  void setup_circle_xz(Scalar & c, const real radius, const real xcent, const real zcent) {  

    for_vijk(c,i,j,k) {
      c[i][j][k] = setup_area(c.xn(i),c.xn(i+1),c.zn(k),c.zn(k+1), xcent, zcent, radius);
      c[i][j][k] /= c.dxc(i)*c.dzc(k);
      if(approx(c[i][j][k],1.0,boil::pico))
        c[i][j][k]=1.0;
    }

    c.exchange_all();
    c.bnd_update();

    return;
  }
}

