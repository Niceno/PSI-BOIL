#include "vof.h"

/******************************************************************************/
void VOF::ev_project(const ScalarInt & pflag, const Matter * fluid,
                     const Scalar & frc, Vector & u) {
/***************************************************************************//**
*  \brief Calculate the velocity correction.
*******************************************************************************/

  Comp m;
  real rho;
  int pfm, pfp;

  m = Comp::u();
  if(ifull)
    for_vmijk(u,m,i,j,k) {
      pfm = abs(pflag[i-1][j][k]);
      pfp = abs(pflag[i  ][j][k]);
      if(pfm<4&&pfp<4&&(pfm+pfp)<6) {
        rho = fluid->rho(m,i,j,k);
        u[m][i][j][k] += (frc[i-1][j][k]-frc[i][j][k]) / rho * time->dt()
                       / (u.dxc(m,i));
      }
    }

  m = Comp::v();
  if(jfull)
    for_vmijk(u,m,i,j,k) {
      pfm = abs(pflag[i][j-1][k]);
      pfp = abs(pflag[i][j  ][k]);
      if(pfm<4&&pfp<4&&(pfm+pfp)<6) {
        rho = fluid->rho(m,i,j,k);
        u[m][i][j][k] += (frc[i][j-1][k]-frc[i][j][k]) / rho * time->dt()
                       / (u.dyc(m,j));
      }
    }

  m = Comp::w();
  if(kfull)
    for_vmijk(u,m,i,j,k) {
      pfm = abs(pflag[i][j][k-1]);
      pfp = abs(pflag[i][j][k  ]);
      if(pfm<4&&pfp<4&&(pfm+pfp)<6) {
        rho = fluid->rho(m,i,j,k);
        u[m][i][j][k] += (frc[i][j][k-1]-frc[i][j][k]) / rho * time->dt()
                       / (u.dzc(m,k));
      }
    }

  u.exchange_all();

  return;
}
