#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::insert_bc_gradt(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief correct gradient of temperature near walls
*******************************************************************************/
  for(int b = 0; b < clr.bc().count(); b++) {

    if(clr.bc().type_decomp(b))
      continue;

    if(clr.bc().type(b) == BndType::wall()) {

      Dir d = clr.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int dir(0), ofx(0), ofy(0), ofz(0);
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          dir = -1;
          ofx = +1;
        } else if (d == Dir::imax()) {
          mcomp = Comp::i();
          dir = +1;
          ofx = -1;
        } else if (d == Dir::jmin()) {
          mcomp = Comp::j();
          dir = -1;
          ofy = +1;
        } else if (d == Dir::jmax()) {
          mcomp = Comp::j();
          dir = +1;
          ofy = -1;
        } else if (d == Dir::kmin()) {
          mcomp = Comp::k();
          dir = -1;
          ofz = +1;
        } else if (d == Dir::kmax()) {
          mcomp = Comp::k();
          dir = +1;
          ofz = -1;
        } else {
          continue;
        }

        for_vijk( clr.bc().at(b), i,j,k ) {
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
          
          /* is there an interface between cell centre and wall? */
          if(Interface(dir,mcomp,ii,jj,kk)) {
            /* wall temperature */
            real tw = tpr[i][j][k];
            /* interface temperature and distance wall-interface */
            real ti;
            real dist;
            /* the temperature gradient for inverse phase is set 
             * extrapolation is turned off in ext_gradt */
            if       (mcomp==Comp::i()) {
              dist = distance_x(ii,jj,kk,dir,ti);
              dist = clr.dxc(ii)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                txv[ii][jj][kk] = (tw-ti)/dist * real(dir);
              } else {
                txl[ii][jj][kk] = (tw-ti)/dist * real(dir);
              }
            } else if(mcomp==Comp::j()) {             
              dist = distance_y(ii,jj,kk,dir,ti);
              dist = clr.dyc(jj)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                tyv[ii][jj][kk] = (tw-ti)/dist * real(dir);
              } else {
                tyl[ii][jj][kk] = (tw-ti)/dist * real(dir);
              }
            } else {
              dist = distance_z(ii,jj,kk,dir,ti);
              dist = clr.dzc(kk)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                tzv[ii][jj][kk] = (tw-ti)/dist * real(dir);
              } else {
                tzl[ii][jj][kk] = (tw-ti)/dist * real(dir);
              }
            }
          }
        }

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */
}
