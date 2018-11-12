#include "levelset.h"

/******************************************************************************/
void LevelSet::insert_bc_gradphic(const Scalar & val) {
/***************************************************************************//**
*  \brief boundary condition for grad(val)
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        real ni,nj,nk,magn;
        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1;
            nx[ii][j][k] = (val[ii+1][j][k]-val[ii  ][j][k])/(      dxe(ii));
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;
            nx[ii][j][k] = (val[ii  ][j][k]-val[ii-1][j][k])/(dxw(ii)      );
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            ny[i][jj][k] = (val[i][jj+1][k]-val[i][jj  ][k])/(      dyn(jj));
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            ny[i][jj][k] = (val[i][jj  ][k]-val[i][jj-1][k])/(dys(jj)      );
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            nz[i][j][kk] = (val[i][j][kk+1]-val[i][j][kk  ])/(      dzt(kk));
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            nz[i][j][kk] = (val[i][j][kk  ]-val[i][j][kk-1])/(dzb(kk)      );
          }
        }
      }
    }

    if( val.bc().type(b) == BndType::symmetry() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        real ni,nj,nk,magn;
        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1; 
            nx[ii][j][k] = (val[ii+1][j][k]-val[ii][j][k])/(2.0*dxw(ii)+dxe(ii));
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1; 
            nx[ii][j][k] = (val[ii][j][k]-val[ii-1][j][k])/(dxw(ii)+2.0*dxe(ii));
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            ny[i][jj][k] = (val[i][jj+1][k]-val[i][jj][k])/(2.0*dys(jj)+dyn(jj));
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            ny[i][jj][k] = (val[i][jj][k]-val[i][jj-1][k])/(dys(j)+2.0*dyn(j));
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            nz[i][j][kk] = (val[i][j][kk+1]-val[i][j][kk])/(2.0*dzb(kk)+dzt(kk));
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            nz[i][j][kk] = (val[i][j][kk]-val[i][j][kk-1])/(dzb(kk)+2.0*dzt(kk));
          }
        }
      }
    }
  }
}
