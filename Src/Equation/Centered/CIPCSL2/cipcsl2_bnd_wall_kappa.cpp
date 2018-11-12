#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::bnd_wall_kappa() {
/***************************************************************************//**
*  \brief Wall boundary condition for curvature.
*  Curvature in wall adjacent cells is used for the interpolation of curvature
*  done in curv_interface().
*  However, the curvature in wall adjacent cells is overwritten by bdcurev.
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<phi.bc().count(); b++ ) {
    if(phi.bc().type_decomp(b)) continue;
    if( phi.bc().type(b) == BndType::wall() ) {
      int iof=0, jof=0, kof=0;
      Dir d      = phi.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        if(d == Dir::imin() || d == Dir::imax()){
          iof=1;
          if(d == Dir::imax()) iof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int ii = i+iof;
            /* x direction */
            if(iof==1){   //  imin
              kappa[ii][j][k] = -(nx[ii+1][j][k]-nx[ii][j][k])/dxe(ii);
            } else {      //  imax
              kappa[ii][j][k] = -(nx[ii][j][k]-nx[ii-1][j][k])/dxw(ii);
            }
            /* y direction */
            kappa[ii][j][k] -= (ny[ii][j+1][k]-ny[ii][j-1][k])/(dys(j)+dyn(j));
            /* z direction */
            kappa[ii][j][k] -= (nz[ii][j][k+1]-nz[ii][j][k-1])/(dzb(k)+dzt(k));
          }
        }

        if(d == Dir::jmin() || d == Dir::jmax()){
          jof=1;
          if(d == Dir::jmax()) jof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int jj=j+jof;
            /* x direction */
            kappa[i][jj][k] = -(nx[i+1][jj][k]-nx[i-1][jj][k])/(dxw(i)+dxe(i));
            /* y direction */
            if(jof==1){  // jmin
              kappa[i][jj][k] -= (ny[i][jj+1][k]-ny[i][jj][k])/dyn(jj);
            } else {     // jmax
              kappa[i][jj][k] -= (ny[i][jj][k]-ny[i][jj-1][k])/dys(jj);
            }
            /* z direction */
            kappa[i][jj][k] -= (nz[i][jj][k+1]-nz[i][jj][k-1])/(dzb(k)+dzt(k));
          }
        }
        if(d == Dir::kmin() || d == Dir::kmax()){
          kof=1;
          if(d == Dir::kmax()) kof=-1;
          for_vijk( phi.bc().at(b), i,j,k ){
            int kk=k+kof;
            /* x direction */
            kappa[i][j][kk] = -(nx[i+1][j][kk]-nx[i-1][j][kk])/(dxw(i)+dxe(i));
            /* y direction */
            kappa[i][j][kk] -= (ny[i][j+1][kk]-ny[i][j-1][kk])/(dys(j)+dyn(j));
            /* z direction */
            if(kof==1){  // kmin
              kappa[i][j][kk] -= (nz[i][j][kk+1]-nz[i][j][kk])/(dzt(kk));
            } else {     // kmax
              kappa[i][j][kk] -= (nz[i][j][kk]-nz[i][j][kk-1])/(dzb(kk));
            }
          }
        }
      }
    }
  }

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    if (dom->ibody().off(i-1,j,k) || dom->ibody().off(i+1,j,k)) {

      /* y direction */
      if(dom->ibody().off(i-1,j,k)){   //  imin
        kappa[i][j][k] = -(nx[i+1][j][k]-nx[i][j][k])/dxe(i);
      } else {      //  imax
        kappa[i][j][k] = -(nx[i][j][k]-nx[i-1][j][k])/dxw(i);
      }
      /* y direction */
      kappa[i][j][k] -= (ny[i][j+1][k]-ny[i][j-1][k])/(dys(j)+dyn(j));
      /* z direction */
      kappa[i][j][k] -= (nz[i][j][k+1]-nz[i][j][k-1])/(dzb(k)+dzt(k));

    } else if (dom->ibody().off(i,j-1,k) || dom->ibody().off(i,j+1,k)) {

      /* x direction */
      kappa[i][j][k] = -(nx[i+1][j][k]-nx[i-1][j][k])/(dxw(i)+dxe(i));
      /* y direction */
      if(dom->ibody().off(i,j-1,k)){  // jmin
        kappa[i][j][k] -= (ny[i][j+1][k]-ny[i][j][k])/dyn(j);
      } else {     // jmax
        kappa[i][j][k] -= (ny[i][j][k]-ny[i][j-1][k])/dys(j);
      }
      /* z direction */
      kappa[i][j][k] -= (nz[i][j][k+1]-nz[i][j][k-1])/(dzb(k)+dzt(k));

    } else if (dom->ibody().off(i,j,k-1) || dom->ibody().off(i,j,k+1)) {

      /* x direction */
      kappa[i][j][k] = -(nx[i+1][j][k]-nx[i-1][j][k]) / (dxw(i)+dxe(i));
      /* y direction */
      kappa[i][j][k] -= (ny[i][j+1][k]-ny[i][j-1][k]) / (dys(j)+dyn(j));
      /* z direction */
      if(dom->ibody().off(i,j,k-1)){  // kmin
        kappa[i][j][k] -= (nz[i][j][k+1]-nz[i][j][k]) / (dzt(k));
      } else {     // kmax
        kappa[i][j][k] -= (nz[i][j][k]-nz[i][j][k-1]) / (dzb(k));
      }

    } 
  }
#endif
}
