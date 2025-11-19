#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::fd_stgd(const Domain & dom, const Body & nfx, const Body & nfy,
                       const Body & nfz) {

#ifdef DEBUG
  std::cout<<"body_fd_stgd:start\n";
#endif

    real fdxw=1.0, fdxe=1.0, 
         fdys=1.0, fdyn=1.0, 
         fdzb=1.0, fdzt=1.0;

    /* "i" direction; dxw & dxe */
    const real dxw = xc - xw;
    const real dxe = xe - xc;
    if( xi < xc && xi > xw ) fdxw = (xc - xi)/dxw;
    if( xi > xc && xi < xe ) fdxe = (xi - xc)/dxe;
    assert( fdxw > 0.0 );
    assert( fdxe > 0.0 );
 
    /* "j" direction; dys & dyn */
    const real dys = yc - ys;
    const real dyn = yn - yc;
    if( yi < yc && yi > ys ) fdys = (yc - yi)/dys;
    if( yi > yc && yi < yn ) fdyn = (yi - yc)/dyn;
    assert( fdys > 0.0 );
    assert( fdyn > 0.0 );

    /* "k" direction; dzb & dzt */
    const real dzb = zc - zb; 
    const real dzt = zt - zc;
    if( zi < zc && zi > zb ) fdzb = (zc - zi)/dzb;
    if( zi > zc && zi < zt ) fdzt = (zi - zc)/dzt;
    assert( fdzb > 0.0 );
    assert( fdzt > 0.0 );

    /* set the values */
    ccell->fdxw(fdxw);
    ccell->fdxe(fdxe);
    ccell->fdys(fdys);
    ccell->fdyn(fdyn);
    ccell->fdzb(fdzb);
    ccell->fdzt(fdzt);

#ifdef DEBUG
  std::cout<<"body_fd_stgd:end\n";
#endif

}
