#include "vof.h"

/******************************************************************************/
void VOF::extend_norm(const Scalar & val) {
/***************************************************************************//**
*  \brief normal vector for cells adjacent to a wall or an immersed boundary
*         obtained using extrapolated volume fractions
*******************************************************************************/

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      /*-------+
      |  Wall  |
      +-------*/

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1;
            norm_mixed(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k], ii,j ,k , val);
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;
            norm_mixed(nx[ii][j][k], ny[ii][j][k], nz[ii][j][k], ii,j ,k , val);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            norm_mixed(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k], i ,jj,k , val);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            norm_mixed(nx[i][jj][k], ny[i][jj][k], nz[i][jj][k], i ,jj,k , val);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            //if(phi[i][j][kk]>boil::pico && phi[i][j][kk]<1. && j==boil::BW)
            //    boil::oout<<"Before "<<i<<" "<<phi[i][j][kk]<<" | "<<nx[i][j][kk]<<" "<<ny[i][j][kk]<<" "<<nz[i][j][kk]<<boil::endl;
            norm_mixed(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk], i ,j ,kk, val);
            //if(phi[i][j][kk]>boil::pico && phi[i][j][kk]<1. && j==boil::BW)
            //    boil::oout<<"After "<<i<<" "<<phi[i][j][kk]<<" | "<<nx[i][j][kk]<<" "<<ny[i][j][kk]<<" "<<nz[i][j][kk]<<boil::endl;
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            norm_mixed(nx[i][j][kk], ny[i][j][kk], nz[i][j][kk], i ,j ,kk, val);

          }
        }
      }
    }

  } /* bcs */

  /***************+
  | immersed body |
  +***************/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    // cell[i][j][k] is wall adjacent cells in fluid domain
    dom->ibody().ijk(cc,&i,&j,&k);

    // west is in solid domain
    if (dom->ibody().off(i-1,j,k)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
      norm_mixed(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, val);
    }
  }

  return;
}
