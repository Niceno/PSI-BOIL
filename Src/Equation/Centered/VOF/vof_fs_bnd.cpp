#include "vof.h"

/******************************************************************************/
void VOF::fs_bnd() {
/***************************************************************************//**
*  \brief Corrects fs at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
******************************************************************************/

  /* tolerance is necessary because of errors */
  /* e.g. 1.0 approx 0.999 */
  real tol = 0.5e-2;

  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( phi.bc().type(b) == BndType::wall() ) {

      Dir d = phi.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int of(0), ofx(0), ofy(0), ofz(0);
#if 0 /* doesn't work without constexpr */
        switch(d) {
          case Dir::imin() : mcomp = Comp::i();
                             of  = +1;
                             ofx = +1;
                             break;  
          case Dir::imax() : mcomp = Comp::i();
                             of  = 0;
                             ofx = -1;
                             break;  
          case Dir::jmin() : mcomp = Comp::j();
                             of  = +1;
                             ofy = +1;
                             break;  
          case Dir::jmax() : mcomp = Comp::j();
                             of  = 0;
                             ofy = -1;
                             break;  
          case Dir::kmin() : mcomp = Comp::k();
                             of  = +1;
                             ofz = +1;
                             break;  
          case Dir::kmax() : mcomp = Comp::k();
                             of  = 0;
                             ofz = -1;
                             break;  
          default : continue;
        }
#else
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          of  = +1;
          ofx = +1;
        } else if (d == Dir::imax()) { 
          mcomp = Comp::i();
          of  = 0;
          ofx = -1;
        } else if (d == Dir::jmin()) { 
          mcomp = Comp::j();
          of  = +1;
          ofy = +1;
        } else if (d == Dir::jmax()) { 
          mcomp = Comp::j();
          of  = 0;
          ofy = -1;
        } else if (d == Dir::kmin()) { 
          mcomp = Comp::k();
          of  = +1;
          ofz = +1;
        } else if (d == Dir::kmax()) { 
          mcomp = Comp::k();
          of  = 0;
          ofz = -1;
        } else {
          continue;
        }
#endif
       
        for_vijk( phi.bc().at(b), i,j,k ) { 
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
          real fsval = fs_val(mcomp,ii,jj,kk);
          bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
          bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
  
          /* interface exists in the dir we are interested in */
          if((flagm && of)||(flagp && !of)) { 
            if(mcomp==Comp::i()) {
              fs[mcomp][i+of][j][k] = phi.xn(ii) + phi.dxc(ii) * fsval;
            } else if(mcomp==Comp::j()) {
              fs[mcomp][i][j+of][k] = phi.yn(jj) + phi.dyc(jj) * fsval;   
            } else {
              fs[mcomp][i][j][k+of] = phi.zn(kk) + phi.dzc(kk) * fsval;
            }
          }
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */
}
