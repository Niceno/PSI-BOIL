#include "momentum.h"

/******************************************************************************/
void Momentum::project_ghost(const Scalar & frc,
                             const Scalar & c,
                             const Scalar & kappa ) {
 
  Comp m;
  real rho;

  real sigma = fluid()->sigma()->value();
  real cint = 0.5;

  /*----------------------------------------------------+
  |  project the velocities on the computational cells  |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |     not the most efficient -> too many divisions    |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |      do not try to introduce some integral form     |
  |       here gradients are constant in each cell      |
  |          and depend on boundary values only         |
  +----------------------------------------------------*/

  if( dom->ibody().nccells() == 0 ) {
    m = Comp::u();
    for_mijk(m,i,j,k) {
        real pw = frc[i-1][j][k];
        real pc = frc[i  ][j][k];
        if ((c[i-1][j][k]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i-1][j][k]+kappa[i][j][k]);
          if (kappa[i-1][j][k]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i-1][j][k] * kappa[i][j][k]
                           / (kappa[i-1][j][k] + kappa[i][j][k]);
          }
          pw += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (pw-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dxc(m,i));
    }

    m = Comp::v();
    for_mijk(m,i,j,k) {
        real ps = frc[i][j-1][k];
        real pc = frc[i][j  ][k];
        if ((c[i][j-1][k]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i][j-1][k]+kappa[i][j][k]);
          if (kappa[i][j-1][k]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i][j-1][k] * kappa[i][j][k]
                           / (kappa[i][j-1][k] + kappa[i][j][k]);
          }
          ps += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (ps-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dyc(m,j));
    }

    m = Comp::w();
    for_mijk(m,i,j,k) {
        real pb = frc[i][j][k-1];
        real pc = frc[i][j][k  ];
        if ((c[i][j][k-1]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i][j][k-1]+kappa[i][j][k]);
          if (kappa[i][j][k-1]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i][j][k-1] * kappa[i][j][k]
                           / (kappa[i][j][k-1] + kappa[i][j][k]);
          }
          pb += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (pb-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dzc(m,k));
    }
  } else {
    m = Comp::u();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i-1,j,k) ) {
        real pw = frc[i-1][j][k];
        real pc = frc[i  ][j][k];
        if ((c[i-1][j][k]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i-1][j][k]+kappa[i][j][k]);
          if (kappa[i-1][j][k]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i-1][j][k] * kappa[i][j][k]
                           / (kappa[i-1][j][k] + kappa[i][j][k]);
          }
          pw += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (pw-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dxc(m,i));
      }
    }
    m = Comp::v();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j-1,k) ) {
        real ps = frc[i][j-1][k];
        real pc = frc[i][j  ][k];
        if ((c[i][j-1][k]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i][j-1][k]+kappa[i][j][k]);
          if (kappa[i][j-1][k]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i][j-1][k] * kappa[i][j][k]
                           / (kappa[i][j-1][k] + kappa[i][j][k]);
          }
          ps += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (ps-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dyc(m,j));
      }
    }
    m = Comp::w();
    for_vmijk(u,m,i,j,k) {
      if( dom->ibody().on_p(i,j,k) && dom->ibody().on_p(i,j,k-1) ) {
        real pb = frc[i][j][k-1];
        real pc = frc[i][j][k  ];
        if ((c[i][j][k-1]-cint)*(c[i][j][k]-cint)<0.0) {
          int iphase = 1;
          if (c[i][j][k]<cint) iphase = -1;
          real kappa_ave = 0.5*(kappa[i][j][k-1]+kappa[i][j][k]);
          if (kappa[i][j][k-1]*kappa[i][j][k]>0.0) {
            kappa_ave = 2.0 * kappa[i][j][k-1] * kappa[i][j][k]
                           / (kappa[i][j][k-1] + kappa[i][j][k]);
          }
          pb += real(iphase)*sigma*kappa_ave;
        }
        u[m][i][j][k] += (pb-pc) * time->dt()
                        / fluid()->rho(m,i,j,k)
                        / (u.dzc(m,k));
      }
    }
  }
  insert_bc();
  u.exchange_all();
}
