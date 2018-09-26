#include "body.h"
#include "../Global/global_approx.h"

/******************************************************************************/
bool Body::cross_seg_y(const int index, 
                       const real xc, const Range<real> y_seg, const real zc,
                       real * y,
                       int  * cut_poly) const {

  /*----------------------------------------+ 
  |  do not browse through all polygons,    |
  |  only throught those stored in polytag  |
  +----------------------------------------*/
  for(int p=0; p<polytags[index].size(); p++) {
    const int c = polytags[index][p];

    /* if plane not parallel to y */
    if( !approx( polys[c].b(), 0.0 ) ) {

      /* y intersection */
      real yi = (-polys[c].a()*xc - polys[c].c()*zc - polys[c].d())
              / polys[c].b();   

      if( y_seg.contains(yi) ) {

        real A0 = area(zc, xc, 
                       polys[c].zn(1), polys[c].xn(1), 
                       polys[c].zn(2), polys[c].xn(2)); 

        real A1 = area(polys[c].zn(0), polys[c].xn(0), 
                       zc, xc, 
                       polys[c].zn(2), polys[c].xn(2));
 
        real A2 = area(polys[c].zn(0), polys[c].xn(0), 
                       polys[c].zn(1), polys[c].xn(1), 
                       zc, xc); 

        bool inside = false;
        if( polys[c].n(1)>-boil::pico && A0>=-boil::pico
                   && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
        if( polys[c].n(1)< boil::pico && A0<= boil::pico
                   && A1<= boil::pico && A2<= boil::pico ) inside=true;

        if(inside) {
          *y        = yi;
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
 '$Id: body_cross_seg_y.cpp,v 1.12 2013/10/02 08:41:35 sato Exp $'/
+-----------------------------------------------------------------------------*/
