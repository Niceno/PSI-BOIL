#include "levelset.h"

/******************************************************************************/
void LevelSet::insert_bc_kappa() {
/***************************************************************************//**
*  \brief Wall boundary condition for grad(val).
*         val: color function of tanh profile
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<phi.bc().count(); b++ ) {

    if(phi.bc().type_decomp(b)) continue;

    if( phi.bc().type(b) == BndType::symmetry() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = phi.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        real ni,nj,nk,magn;
        if(d == Dir::imin() || d == Dir::imax()){
          iof=1;
          if(d == Dir::imax()) iof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int ii = i+iof;
            real nw = nx[ii-1][j][k];
            real ne = nx[ii+1][j][k];
            real ns = ny[ii][j-1][k];
            real nn = ny[ii][j+1][k];
            real nb = nz[ii][j][k-1];
            real nt = nz[ii][j][k+1];
            real dxx;
            if(iof==1){
              dxx=2.0*phi.dxw(ii)+phi.dxe(ii);
            } else {
              dxx=phi.dxw(ii)+2.0*phi.dxe(ii);
            }
            real dyy=phi.dys(j)+phi.dyn(j);
            if(j==sj()) dyy = 2.0*phi.dys(j)+phi.dyn(j);
            if(j==ej()) dyy = phi.dys(j)+2.0*phi.dyn(j);
            real dzz=phi.dzb(k)+phi.dzt(k);
            if(k==sk()) dzz = 2.0*phi.dzb(k)+phi.dzt(k);
            if(k==ek()) dzz = phi.dzb(k)+2.0*phi.dzt(k);
            kappa[ii][j][k]=-((ne-nw)/dxx
                             +(nn-ns)/dyy
                             +(nt-nb)/dzz);
          }
        }

        if(d == Dir::jmin() || d == Dir::jmax()){
          jof=1;
          if(d == Dir::jmax()) jof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int jj=j+jof;
            real dxx=phi.dxw(i)+phi.dxe(i);
            if(i==si()) dxx = 2.0*phi.dxw(i)+phi.dxe(i);
            if(i==ei()) dxx = phi.dxw(i)+2.0*phi.dxe(i);
            real dyy;
            if(jof==1){
              dyy = 2.0*phi.dys(jj)+phi.dyn(jj);
            } else {
              dyy = phi.dys(jj)+2.0*phi.dyn(jj);
            }
            real dzz=phi.dzb(k)+phi.dzt(k);
            if(k==sk()) dzz = 2.0*phi.dzb(k)+phi.dzt(k);
            if(k==ek()) dzz = phi.dzb(k)+2.0*phi.dzt(k);
            real nw = nx[i-1][jj][k];
            real ne = nx[i+1][jj][k];
            real ns = ny[i][jj-1][k];
            real nn = ny[i][jj+1][k];
            real nb = nz[i][jj][k-1];
            real nt = nz[i][jj][k+1];
            kappa[i][jj][k]=-((ne-nw)/dxx
                             +(nn-ns)/dyy
                             +(nt-nb)/dzz);
          }
        }
        if(d == Dir::kmin() || d == Dir::kmax()){
          kof=1;
          if(d == Dir::kmax()) kof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int kk=k+kof;
            real dxx=phi.dxw(i)+phi.dxe(i);
            if(i==si()) dxx = 2.0*phi.dxw(i)+phi.dxe(i);
            if(i==ei()) dxx = phi.dxw(i)+2.0*phi.dxe(i);
            real dyy=phi.dys(j)+phi.dyn(j);
            if(j==sj()) dyy = 2.0*phi.dys(j)+phi.dyn(j);
            if(j==ej()) dyy = phi.dys(j)+2.0*phi.dyn(j);
            real dzz;
            if(kof==1){
              dzz=2.0*phi.dzb(kk)+phi.dzt(kk);
            } else {
              dzz=phi.dzb(kk)+2.0*phi.dzt(kk);
            }
            real nw = nx[i-1][j][kk];
            real ne = nx[i+1][j][kk];
            real ns = ny[i][j-1][kk];
            real nn = ny[i][j+1][kk];
            real nb = nz[i][j][kk-1];
            real nt = nz[i][j][kk+1];
            kappa[i][j][kk]=-((ne-nw)/dxx
                             +(nn-ns)/dyy
                             +(nt-nb)/dzz);
          }
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: levelset_insert_bc_kappa.cpp,v 1.2 2012/09/13 08:42:27 niceno Exp $'/
+-----------------------------------------------------------------------------*/
