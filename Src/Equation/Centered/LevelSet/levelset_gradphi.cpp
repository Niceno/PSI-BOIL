#include "levelset.h"

/******************************************************************************/
void LevelSet::gradphi() {
/***************************************************************************//**
*  \brief Calculate grad(csa)/ at node point.
*         Reference: J.U.Brackbill, et al.,"A Continum method for modeling
*                    surface tension",J.Comp.phys.,Vol.100,pp.335-354,1992
*                    Equation (41)
*         Resluts: nx, ny, nz
*******************************************************************************/

  real ni,nj,nk,magn;

  /* node base */
  for(int i=phi.si(); i<=phi.ei()+1; i++) {
  for(int j=phi.sj(); j<=phi.ej()+1; j++) {
  for(int k=phi.sk(); k<=phi.ek()+1; k++) {
        real dx=dxw(i);
        real dy=dys(j);
        real dz=dzb(k);
        nx[i][j][k] = 0.25*( (phi[i][j  ][k  ]-phi[i-1][j  ][k  ])/dx
                            +(phi[i][j-1][k  ]-phi[i-1][j-1][k  ])/dx
                            +(phi[i][j  ][k-1]-phi[i-1][j  ][k-1])/dx
                            +(phi[i][j-1][k-1]-phi[i-1][j-1][k-1])/dx);
        ny[i][j][k] = 0.25*( (phi[i  ][j][k  ]-phi[i  ][j-1][k  ])/dy
                            +(phi[i-1][j][k  ]-phi[i-1][j-1][k  ])/dy
                            +(phi[i  ][j][k-1]-phi[i  ][j-1][k-1])/dy
                            +(phi[i-1][j][k-1]-phi[i-1][j-1][k-1])/dy);
        nz[i][j][k] = 0.25*( (phi[i  ][j  ][k]-phi[i  ][j  ][k-1])/dz
                            +(phi[i  ][j-1][k]-phi[i  ][j-1][k-1])/dz
                            +(phi[i-1][j  ][k]-phi[i-1][j  ][k-1])/dz
                            +(phi[i-1][j-1][k]-phi[i-1][j-1][k-1])/dz);
  } } }
  
  insert_bc_gradphi(phi);

  for(int i=phi.si(); i<=phi.ei()+1; i++) {
  for(int j=phi.sj(); j<=phi.ej()+1; j++) {
  for(int k=phi.sk(); k<=phi.ek()+1; k++) {
      ni = nx[i][j][k];
      nj = ny[i][j][k];
      nk = nz[i][j][k];
      real magn = sqrt(ni*ni + nj*nj + nk*nk) + epsnorm;
      nx[i][j][k] = ni/magn;
      ny[i][j][k] = nj/magn;
      nz[i][j][k] = nk/magn;
  } } }

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: levelset_gradphi.cpp,v 1.2 2012/09/13 08:42:26 niceno Exp $'/
+-----------------------------------------------------------------------------*/
