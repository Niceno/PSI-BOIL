#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_gradphic(const Scalar & val) {
/***************************************************************************//**
*  \brief normal vector for the adjacent cells next wall and immersed boundary
*   Note: values are set inside the domain, not in dummy cells!
*******************************************************************************/

  int i,j,k;

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
            if (i<=si()-2) continue;
            if (i>=ei()+2) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            int ii=i+1;
            nx[ii][j][k] = (val[ii+1][j][k]-val[ii  ][j][k])/(      dxe(ii));
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            if (i<=si()-2) continue;
            if (i>=ei()+2) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            int ii=i-1;
            nx[ii][j][k] = (val[ii  ][j][k]-val[ii-1][j][k])/(dxw(ii)      );
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-2) continue;
            if (j>=ej()+2) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            int jj=j+1;
            ny[i][jj][k] = (val[i][jj+1][k]-val[i][jj  ][k])/(      dyn(jj));
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-2) continue;
            if (j>=ej()+2) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            int jj=j-1;
            ny[i][jj][k] = (val[i][jj  ][k]-val[i][jj-1][k])/(dys(jj)      );
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-2) continue;
            if (k>=ek()+2) continue;
            int kk=k+1;
            nz[i][j][kk] = (val[i][j][kk+1]-val[i][j][kk  ])/(      dzt(kk));
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-2) continue;
            if (k>=ek()+2) continue;
            int kk=k-1;
            nz[i][j][kk] = (val[i][j][kk  ]-val[i][j][kk-1])/(dzb(kk)      );
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
      nx[i][j][k]=(val[i+1][j][k]-val[i][j][k])/dxe(i);
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
      nx[i][j][k]=(val[i][j][k]-val[i-1][j][k])/dxw(i);
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
      ny[i][j][k]=(val[i][j+1][k]-val[i][j][k])/dyn(j);
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
      ny[i][j][k]=(val[i][j][k]-val[i][j-1][k])/dys(j);
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
     nz[i][j][k]=(val[i][j][k+1]-val[i][j][k])/dzt(k);
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
      nz[i][j][k]=(val[i][j][k]-val[i][j][k-1])/dzb(k);
    }
  }
#endif
}
