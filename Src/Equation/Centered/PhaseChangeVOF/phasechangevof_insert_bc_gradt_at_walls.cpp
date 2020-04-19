#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::insert_bc_gradt_at_walls(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief correct gradient of temperature at walls
*         N.B. fluid values are also updated but will be overwritten during
*         extrapolation
*******************************************************************************/
  for(int b = 0; b < tpr.bc().count(); b++) {

    if(tpr.bc().type_decomp(b))
      continue;

    if(tpr.bc().type(b) == BndType::dirichlet()) {

      Dir d = clr.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int ofx(0), ofy(0), ofz(0);
        Sign sig;
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          sig = Sign::neg();
          ofx = +1;
        } else if (d == Dir::imax()) {
          mcomp = Comp::i();
          sig = Sign::pos();
          ofx = -1;
        } else if (d == Dir::jmin()) {
          mcomp = Comp::j();
          sig = Sign::neg();
          ofy = +1;
        } else if (d == Dir::jmax()) {
          mcomp = Comp::j();
          sig = Sign::pos();
          ofy = -1;
        } else if (d == Dir::kmin()) {
          mcomp = Comp::k();
          sig = Sign::neg();
          ofz = +1;
        } else if (d == Dir::kmax()) {
          mcomp = Comp::k();
          sig = Sign::pos();
          ofz = -1;
        } else {
          continue;
        }

        for_vijk( tpr.bc().at(b), i,j,k ) {
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
          
          /* is there an interface between cell centre and wall? */
          if(Interface(sig,mcomp,ii,jj,kk)) {

            /* wall temperature */
            real tw = tpr[i][j][k];
            /* interface temperature and distance wall-interface */
            real ti;
            real dist;
            /* the temperature gradient for inverse phase is set 
             * extrapolation overwrites them though */
            if       (mcomp==Comp::i()) {
              dist = distance_x(ii,jj,kk,sig,ti);
              dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                txv[ii][jj][kk] = (tw-ti)/dist * real(sig);
                txv[i ][j ][k ] = (tw-ti)/dist * real(sig);
              } else {
                txl[ii][jj][kk] = (tw-ti)/dist * real(sig);
                txl[i ][j ][k ] = (tw-ti)/dist * real(sig);
              }
            } else if(mcomp==Comp::j()) {             
              dist = distance_y(ii,jj,kk,sig,ti);
              dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                tyv[ii][jj][kk] = (tw-ti)/dist * real(sig);
                tyv[i ][j ][k ] = (tw-ti)/dist * real(sig);
              } else {
                tyl[ii][jj][kk] = (tw-ti)/dist * real(sig);
                tyl[i ][j ][k ] = (tw-ti)/dist * real(sig);
              }
            } else {
              dist = distance_z(ii,jj,kk,sig,ti);
              dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                tzv[ii][jj][kk] = (tw-ti)/dist * real(sig);
                tzv[i ][j ][k ] = (tw-ti)/dist * real(sig);
              } else {
                tzl[ii][jj][kk] = (tw-ti)/dist * real(sig);
                tzl[i ][j ][k ] = (tw-ti)/dist * real(sig);
              }
            }

          } /* there is interface */
        } /* ijk */

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  return;
}
