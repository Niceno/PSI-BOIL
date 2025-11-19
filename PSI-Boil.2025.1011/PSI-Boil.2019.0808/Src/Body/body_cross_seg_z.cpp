#include "body.h"
#include "../Global/global_approx.h"

/******************************************************************************/
bool Body::cross_seg_z(const int index, 
                       const real xc, const real yc, const Range<real> z_seg, 
                       real * z,
                       int  * cut_poly) const {

  /*----------------------------------------+ 
  |  do not browse through all polygons,    |
  |  only throught those stored in polytag  |
  +----------------------------------------*/
  for(int p=0; p<polytags[index].size(); p++) {
    const int c = polytags[index][p];

    /* if plane not parallel to z */
    if( !approx( polys[c].c(), 0.0 ) ) {

      /* z intersection */
      real zi = (-polys[c].a()*xc - polys[c].b()*yc - polys[c].d())
              / polys[c].c();   

      if( z_seg.contains(zi) ) {

        real A0 = area(xc, yc, 
                       polys[c].xn(1), polys[c].yn(1), 
                       polys[c].xn(2), polys[c].yn(2)); 

        real A1 = area(polys[c].xn(0), polys[c].yn(0), 
                       xc, yc, 
                       polys[c].xn(2), polys[c].yn(2));
 
        real A2 = area(polys[c].xn(0), polys[c].yn(0), 
                       polys[c].xn(1), polys[c].yn(1), 
                       xc, yc); 

        bool inside = false;
        if( polys[c].n(2)>-boil::pico && A0>=-boil::pico
                   && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
        if( polys[c].n(2)< boil::pico && A0<= boil::pico
                   && A1<= boil::pico && A2<= boil::pico ) inside=true;


        if(inside) {
          *z        = zi;
          *cut_poly = c;

          return true;
        }
      }
    }
  }  

  *cut_poly = -1;
  return false;
}
