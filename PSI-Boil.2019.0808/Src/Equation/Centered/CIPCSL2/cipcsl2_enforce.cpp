#include "cipcsl2.h"
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::enforce ( const int i, const int j, const int k, const real xx) {
/***************************************************************************//**
*  \brief enforce a value to: phi, clr
*******************************************************************************/

  phi[i][j][k]=xx;
  clr[i][j][k]=xx;

  /* f-node */
  for(int ii=0; ii<=1; ii++){
  for(int jj=0; jj<=1; jj++){
  for(int kk=0; kk<=1; kk++){
    scheme.f[i+ii][j+jj][k+kk]=xx;
  }}}

  /* line density sigx */
  real dx=phi.dxc(i);
  for(int jj=0; jj<=1; jj++){
  for(int kk=0; kk<=1; kk++){
    scheme.sigx[i][j+jj][k+kk]=xx*dx;
  }}

  /* line density sigy */
  real dy=phi.dyc(j);
  for(int ii=0; ii<=1; ii++){
  for(int kk=0; kk<=1; kk++){
    scheme.sigx[i+ii][j][k+kk]=xx*dy;
  }}

  /* line density sigz */
  real dz=phi.dzc(k);
  for(int ii=0; ii<=1; ii++){
  for(int jj=0; jj<=1; jj++){
    scheme.sigz[i+ii][j+jj][k]=xx*dz;
  }}

  /* face density sxyz[i] */
  for(int ii=0; ii<=1; ii++){
    sxyz[Comp::i()][i+ii][j][k]=xx*dz*dy;
  }

  /* face density sxyz[j] */
  for(int jj=0; jj<=1; jj++){
    sxyz[Comp::j()][i][j+jj][k]=xx*dz*dx;
  }

  /* face density sxyz[k] */
  for(int kk=0; kk<=1; kk++){
    sxyz[Comp::k()][i][j][k+kk]=xx*dx*dy;
  }

  return;
}
