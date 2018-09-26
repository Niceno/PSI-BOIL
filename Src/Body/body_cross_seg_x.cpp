#include "body.h"
#include "../Global/global_approx.h"

/******************************************************************************/
bool Body::cross_seg_x(const int index, 
                       const Range<real> x_seg, const real yc, const real zc,
                       real * x,
                       int  * cut_poly) const {

  /*----------------------------------------+ 
  |  do not browse through all polygons,    |
  |  only throught those stored in polytag  |
  +----------------------------------------*/
  for(int p=0; p<polytags[index].size(); p++) {
    const int c = polytags[index][p];

    /* if plane not parallel to x */
    if( !approx( polys[c].a(), 0.0 ) ) {

      /* x intersection */
      real xi = (-polys[c].b()*yc - polys[c].c()*zc - polys[c].d())
              / polys[c].a();

      if( x_seg.contains(xi) ) {
 
        real A0 = area(yc, zc, 
                       polys[c].yn(1), polys[c].zn(1), 
                       polys[c].yn(2), polys[c].zn(2)); 

        real A1 = area(polys[c].yn(0), polys[c].zn(0), 
                       yc, zc, 
                       polys[c].yn(2), polys[c].zn(2));
 
        real A2 = area(polys[c].yn(0), polys[c].zn(0), 
                       polys[c].yn(1), polys[c].zn(1), 
                       yc, zc); 

        bool inside = false;
        if( polys[c].n(0)>-boil::pico && A0>=-boil::pico
                   && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
        if( polys[c].n(0)< boil::pico && A0<= boil::pico
                   && A1<= boil::pico && A2<= boil::pico ) inside=true;

        if(inside) {
          *x        = xi;
          *cut_poly = c;

          return true;
        }
      }
    }
  }  

  *cut_poly = -1;
  return false;
}

/*-----------------------------------------------------------------------------+
 '$Id: body_cross_seg_x.cpp,v 1.12 2013/10/02 08:41:35 sato Exp $'/
+-----------------------------------------------------------------------------*/
