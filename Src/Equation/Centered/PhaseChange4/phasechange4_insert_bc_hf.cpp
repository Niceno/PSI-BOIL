#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::insert_bc_hf(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate heat flux in subgrid films at boundaries with dirichlet bc.
*         N.B. fluid values are also updated but will be overwritten during
*         extrapolation, since extrapolation also 'looks' into buffer cells.
*******************************************************************************/
  for(int b = 0; b < cht.tmp().bc().count(); b++) {

    if(cht.tmp().bc().type_decomp(b))
      continue;

    if(cht.tmp().bc().type(b) == BndType::dirichlet()) {

      Dir d = phi.bc().direction(b);
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

        for_vijk( cht.tmp().bc().at(b), i,j,k ) {
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
          
          /* is there an interface between cell centre and wall? */
          if(cht.interface(sig,mcomp,ii,jj,kk)) {

            /* wall temperature */
            real tw = cht.tmp()[i][j][k];
            /* interface temperature and distance wall-interface */
            real ti;
            real dist;
            /* subgrid thermal conductivity */
            real lmb = cht.lambda_inv(ii,jj,kk,diff_eddy);
            /* the temperature gradient for inverse phase is set 
             * extrapolation overwrites them though */
            if       (mcomp==Comp::i()) {
              //dist = distance_int_x(sig,ii,jj,kk,ti);
              //dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              dist = cht.distance_int_x(-sig,i,j,k,ti);
              if(cht.topo->above_interface(ii,jj,kk)) {
                txv[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                txv[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              } else {
                txl[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                txl[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              }
            } else if(mcomp==Comp::j()) {             
              //dist = distance_int_y(sig,ii,jj,kk,ti);
              //dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              dist = cht.distance_int_y(-sig,i,j,k,ti);
              if(cht.topo->above_interface(ii,jj,kk)) {
                tyv[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                tyv[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              } else {
                tyl[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                tyl[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              }
            } else {
              //dist = distance_int_z(sig,ii,jj,kk,ti);
              //dist = distance_center(sig,mcomp,ii,jj,kk) - dist;
              dist = cht.distance_int_z(-sig,i,j,k,ti);
              if(cht.topo->above_interface(ii,jj,kk)) {
                tzv[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                tzv[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              } else {
                tzl[ii][jj][kk] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
                tzl[i ][j ][k ] = (tw-ti)/(dist/lmb+cht.wall_resistance(ii,jj,kk))*real(sig);
              }
            }

          } /* there is interface */
        } /* ijk */

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  return;
}
