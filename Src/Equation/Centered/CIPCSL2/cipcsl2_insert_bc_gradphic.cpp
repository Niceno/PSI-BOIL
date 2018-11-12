#include "cipcsl2.h"
#define ORDER_1st

/******************************************************************************/
void CIPCSL2::insert_bc_gradphic(const Scalar & val) {
/***************************************************************************//**
*  \brief normal vector for the adjacent cells next wall, symmetric and
*         immersed boundary
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
#ifdef ORDER_1st
            nx[ii][j][k] = (val[ii+1][j][k]-val[ii  ][j][k])/(      dxe(ii));
#else
            real a = dxe(ii);
            real b = dxe(ii+1);
            nz[ii][j][k] = -b*(2.0*a+b)*val[ii  ][j][k]
                         + (a+b)*(a+b) *val[ii+1][j][k]
                         - a*a         *val[ii+2][j][k];
            nz[ii][j][k] /= a*b*(a+b);
#endif
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;
#ifdef ORDER_1st
            nx[ii][j][k] = (val[ii  ][j][k]-val[ii-1][j][k])/(dxw(ii)      );
#else
            real a = dxw(ii);
            real b = dxw(ii-1);
            nz[ii][j][k] = -b*(2.0*a+b)*val[ii  ][j][k]
                         + (a+b)*(a+b) *val[ii-1][j][k]
                         - a*a         *val[ii-2][j][k];
            nz[ii][j][k] /= a*b*(a+b);
#endif
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
#ifdef ORDER_1st
            ny[i][jj][k] = (val[i][jj+1][k]-val[i][jj  ][k])/(      dyn(jj));
#else
            real a = dyn(jj);
            real b = dyn(jj+1);
            nz[i][jj][k] = -b*(2.0*a+b)*val[i][jj  ][k]
                         + (a+b)*(a+b) *val[i][jj+1][k]
                         - a*a         *val[i][jj+2][k];
            nz[i][jj][k] /= a*b*(a+b);
#endif
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
#ifdef ORDER_1st
            ny[i][jj][k] = (val[i][jj  ][k]-val[i][jj-1][k])/(dys(jj)      );
#else
            real a = dys(jj);
            real b = dys(jj-1);
            nz[i][jj][k] = -b*(2.0*a+b)*val[i][jj  ][k]
                         + (a+b)*(a+b) *val[i][jj-1][k]
                         - a*a         *val[i][jj-2][k];
            nz[i][jj][k] /= a*b*(a+b);
#endif
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
#ifdef ORDER_1st
            nz[i][j][kk] = (val[i][j][kk+1]-val[i][j][kk  ])/(      dzt(kk));
#else
            real a = dzt(kk);
            real b = dzt(kk+1);
            nz[i][j][kk] = -b*(2.0*a+b)*val[i][j][kk]
                         + (a+b)*(a+b) *val[i][j][kk+1]
                         - a*a         *val[i][j][kk+2];
            nz[i][j][kk] /= a*b*(a+b);
#endif
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
#ifdef ORDER_1st
            nz[i][j][kk] = (val[i][j][kk  ]-val[i][j][kk-1])/(dzb(kk)      );
#else
            real a = dzb(kk);
            real b = dzb(kk-1);
            nz[i][j][kk] = -b*(2.0*a+b)*val[i][j][kk]
                         + (a+b)*(a+b) *val[i][j][kk-1]
                         - a*a         *val[i][j][kk-2];
            nz[i][j][kk] /= a*b*(a+b);
#endif
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

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    // west
    if (dom->ibody().off(i-1,j,k)) {
#ifdef ORDER_1st
      nx[i][j][k]=(val[i+1][j][k]-val[i][j][k])/dxe(i);
#else
      if (i+2>ei()+1) {
        nx[i][j][k]=(val[i+1][j][k]-val[i][j][k])/dxe(i);
      } else {
        real a = dxe(i);
        real b = dxe(i+1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i  ][j][k]
                    + (a+b)*(a+b) *val[i+1][j][k]
                    - a*a         *val[i+2][j][k];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
#ifdef ORDER_1st
      nx[i][j][k]=(val[i][j][k]-val[i-1][j][k])/dxw(i);
#else
      if (i-2<si()-1) {
        nx[i][j][k]=(val[i][j][k]-val[i-1][j][k])/dxw(i);
      } else {
        real a = dxw(i);
        real b = dxw(i-1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i  ][j][k]
                    + (a+b)*(a+b) *val[i-1][j][k]
                    - a*a         *val[i-2][j][k];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
#ifdef ORDER_1st
      ny[i][j][k]=(val[i][j+1][k]-val[i][j][k])/dyn(j);
#else
      if (j+2>ej()+1) {
        ny[i][j][k]=(val[i][j+1][k]-val[i][j][k])/dyn(j);
      } else {
        real a = dyn(j);
        real b = dyn(j+1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i][j  ][k]
                    + (a+b)*(a+b) *val[i][j+1][k]
                    - a*a         *val[i][j+2][k];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
#ifdef ORDER_1st
      ny[i][j][k]=(val[i][j][k]-val[i][j-1][k])/dys(j);
#else
      if (j-2<sj()-1) {
        ny[i][j][k]=(val[i][j][k]-val[i][j-1][k])/dys(j);
      } else {
        real a = dys(j);
        real b = dys(j-1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i][j  ][k]
                    + (a+b)*(a+b) *val[i][j-1][k]
                    - a*a         *val[i][j-2][k];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
#ifdef ORDER_1st
     nz[i][j][k]=(val[i][j][k+1]-val[i][j][k])/dzt(k);
#else
      if (k+2>ek()+1) {                                   // crude code
        nz[i][j][k]=(val[i][j][k+1]-val[i][j][k])/dzt(k);
      } else {
        real a = dzt(k);
        real b = dzt(k+1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i][j][k]
                    + (a+b)*(a+b) *val[i][j][k+1]
                    - a*a         *val[i][j][k+2];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
#ifdef ORDER_1st
      nz[i][j][k]=(val[i][j][k]-val[i][j][k-1])/dzb(k);
#else
      if (k-2<sk()-1) {                                   // crude code
        nz[i][j][k]=(val[i][j][k]-val[i][j][k-1])/dzb(k);
      } else {
        real a = dzb(k);
        real b = dzb(k-1);
        nz[i][j][k] = -b*(2.0*a+b)*val[i][j][k]
                     + (a+b)*(a+b)*val[i][j][k-1]
                     - a*a        *val[i][j][k-2];
        nz[i][j][k] /= a*b*(a+b);
      }
#endif
    }
  }
#endif
}
