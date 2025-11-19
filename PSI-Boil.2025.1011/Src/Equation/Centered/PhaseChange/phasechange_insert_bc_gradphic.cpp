#include "phasechange.h"

/******************************************************************************/
void PhaseChange::insert_bc_gradphic(const Scalar & val) {
/***************************************************************************//**
*  \brief boundary condition for the adjacent cell next to boundary
*******************************************************************************/

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
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

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1; 
            nx[ii][j][k] = (val[ii+1][j][k]-val[ii][j][k])/(dxw(ii)+dxe(ii));
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1; 
            nx[ii][j][k] = (val[ii][j][k]-val[ii-1][j][k])/(dxw(ii)+dxe(ii));
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            ny[i][jj][k] = (val[i][jj+1][k]-val[i][jj][k])/(dys(jj)+dyn(jj));
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            ny[i][jj][k] = (val[i][jj][k]-val[i][jj-1][k])/(dys(j)+dyn(j));
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            nz[i][j][kk] = (val[i][j][kk+1]-val[i][j][kk])/(dzb(kk)+dzt(kk));
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            nz[i][j][kk] = (val[i][j][kk]-val[i][j][kk-1])/(dzb(kk)+dzt(kk));
          }
        }
      }
    }
  }

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    // w
    if(dom->ibody().off(i-1,j,k)){
      nx[i][j][k]=(val[i+1][j][k]-val[i][j][k])/dxe(i);
    }
    // e
    if(dom->ibody().off(i+1,j,k)){
      nx[i][j][k]=(val[i][j][k]-val[i-1][j][k])/dxw(i);
    }
    // s
    if(dom->ibody().off(i,j-1,k)){
      ny[i][j][k]=(val[i][j+1][k]-val[i][j][k])/dyn(j);
    }
    // n
    if(dom->ibody().off(i,j+1,k)){
      ny[i][j][k]=(val[i][j][k]-val[i][j-1][k])/dys(j);
    }
    // b
    if(dom->ibody().off(i,j,k-1)){
      nz[i][j][k]=(val[i][j][k+1]-val[i][j][k])/dzt(k);
    }
    // t
    if(dom->ibody().off(i,j,k+1)){
      nz[i][j][k]=(val[i][j][k]-val[i][j][k-1])/dzb(k);
    }
  }
#endif
}
