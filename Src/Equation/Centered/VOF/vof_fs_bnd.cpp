#include "vof.h"

/******************************************************************************/
void VOF::fs_bnd() {
/***************************************************************************//**
*  \brief Corrects fs at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
*         IMPORTANT: does not work when immersed boundaries do not correspond
*                    to cell boundaries!!!
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

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in solid domain */
    if (dom->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i  ][j][k] = phi.xn(i) + phi.dxc(i) * fsval;
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i+1][j][k] = phi.xn(i) + phi.dxc(i) * fsval;
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j  ][k] = phi.yn(j) + phi.dyc(j) * fsval;
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i][j+1][k] = phi.yn(j) + phi.dyc(j) * fsval;
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j][k  ] = phi.zn(k) + phi.dzc(k) * fsval;
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i][j][k+1] = phi.zn(k) + phi.dzc(k) * fsval;
    }
  }

}
