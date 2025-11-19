#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::bdphiface(const Vector & vec
                      , const Comp & m
                      , const Scalar & sca) {
/***************************************************************************//**
*  \brief Set boundary condition for face variables.
*    exchange() or exchange_all() should take account of periodic condition.
*    1st: bdphiface(vec,m,sca)
*    2nd: vec.exchange_all()
*******************************************************************************/

  /* i-min or i-max */
  int iof,iof2;
  if(m==Comp::i()){
    for( int b=0; b<sca.bc().count(); b++ ) {

      if( sca.bc().type_decomp(b) ) continue;

      /* i-min or i-max */
      if ( phi.bc().direction(b)==Dir::imin()
        || phi.bc().direction(b)==Dir::imax() ){
        if (phi.bc().direction(b)==Dir::imin()){
          iof=1; iof2=1;
        }
        if (phi.bc().direction(b)==Dir::imax()){
          iof=-1; iof2=0;
        }
        /* inlet or dirichlet */
        if( phi.bc().type(b) == BndType::inlet()
          ||phi.bc().type(b) == BndType::dirichlet()
          ||phi.bc().type(b) == BndType::insert() ){
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-2) continue;
            if (i>=ei()+2) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            real dy=phi.dyc(j);
            real dz=phi.dzc(k);
            real stmp = min(1.0,max(0.0,sca[i][j][k]));
            vec[m][i+iof2][j][k] = dy*dz*stmp;
          }
        }
        /* others */
        if( phi.bc().type(b) == BndType::neumann()
          ||phi.bc().type(b) == BndType::symmetry()
          ||phi.bc().type(b) == BndType::wall()
          ||phi.bc().type(b) == BndType::outlet() ) {
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-2) continue;
            if (i>=ei()+2) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            real dy=phi.dyc(j);
            real dz=phi.dzc(k);
            real stmp = min(1.0,max(0.0,sca[i+iof][j][k]));
            vec[m][i+iof2][j][k] = dy*dz*stmp;
          }
        }
      }
    }
  }

  /* j-min or j-max */
  int jof,jof2;
  if(m==Comp::j()){
    for( int b=0; b<sca.bc().count(); b++ ) {

      if( sca.bc().type_decomp(b) ) continue;

      /* j-min or j-max */
      if ( phi.bc().direction(b)==Dir::jmin()
        || phi.bc().direction(b)==Dir::jmax() ){
        if (phi.bc().direction(b)==Dir::jmin()){
          jof=1; jof2=1;
        }
        if (phi.bc().direction(b)==Dir::jmax()){
          jof=-1; jof2=0;
        }
        /* inlet or dirichlet */
        if( phi.bc().type(b) == BndType::inlet()
          ||phi.bc().type(b) == BndType::dirichlet()
          ||phi.bc().type(b) == BndType::insert() ){
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-2) continue;
            if (j>=ej()+2) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            real dx=phi.dxc(i);
            real dz=phi.dzc(k);
            real stmp = min(1.0,max(0.0,sca[i][j][k]));
            vec[m][i][j+jof2][k] = dx*dz*stmp;
          }
        }
        /* others */
        if( phi.bc().type(b) == BndType::neumann()
          ||phi.bc().type(b) == BndType::symmetry()
          ||phi.bc().type(b) == BndType::wall()
          ||phi.bc().type(b) == BndType::outlet() ) {
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-2) continue;
            if (j>=ej()+2) continue;
            if (k<=sk()-1) continue;
            if (k>=ek()+1) continue;
            real dx=phi.dxc(i);
            real dz=phi.dzc(k);
            real stmp = min(1.0,max(0.0,sca[i][j+jof][k]));
            vec[m][i][j+jof2][k] = dx*dz*stmp;
          }
        }
      }
    }
  }

  /* k-min or k-max */
  int kof,kof2;
  if(m==Comp::k()){
    for( int b=0; b<sca.bc().count(); b++ ) {

      if( sca.bc().type_decomp(b) ) continue;

      if ( phi.bc().direction(b)==Dir::kmin()
        || phi.bc().direction(b)==Dir::kmax() ){
        if (phi.bc().direction(b)==Dir::kmin()){
          kof=1; kof2=1;
        }
        if (phi.bc().direction(b)==Dir::kmax()){
          kof=-1; kof2=0;
        }
        /* inlet or dirichlet */
        if( phi.bc().type(b) == BndType::inlet()
          ||phi.bc().type(b) == BndType::dirichlet()
          ||phi.bc().type(b) == BndType::insert() ){
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-2) continue;
            if (k>=ek()+2) continue;
            real dx=phi.dxc(i);
            real dy=phi.dyc(j);
            real stmp = min(1.0,max(0.0,sca[i][j][k]));
            vec[m][i][j][k+kof2] = dx*dy*stmp;
          }
        }
        /* others */
        if( phi.bc().type(b) == BndType::neumann()
          ||phi.bc().type(b) == BndType::symmetry()
          ||phi.bc().type(b) == BndType::wall()
          ||phi.bc().type(b) == BndType::outlet() ) {
          for_vijk( phi.bc().at(b), i,j,k ) {
            if (i<=si()-1) continue;
            if (i>=ei()+1) continue;
            if (j<=sj()-1) continue;
            if (j>=ej()+1) continue;
            if (k<=sk()-2) continue;
            if (k>=ek()+2) continue;
            real dx=phi.dxc(i);
            real dy=phi.dyc(j);
            real stmp = min(1.0,max(0.0,sca[i][j][k+kof]));
            vec[m][i][j][k+kof2] = dx*dy*stmp;
          }
        }
      }
    }
  }
}
