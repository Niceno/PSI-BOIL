#include "cipcsl2.h"
#include <iomanip>
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::init() {
/***************************************************************************//**
*  \brief  Initialize variables used for CIPCSL2
*******************************************************************************/

  for_aijk(i,j,k){
    phi[i][j][k]=min(1.0,max(0.0,phi[i][j][k]));
  }
 
  bdcond(phi);
  phi.exchange_all();

  /* initial condition of color-function:clr */
  for_aijk(i,j,k){
    clr[i][j][k]=phi[i][j][k];
  }

  for_ijk(i,j,k){
    if(fabs(clr[i][j][k]-phisurf)<eps_clr) {
      clr[i][j][k]=phisurf+copysign(1.0,phi[i][j][k]-phisurf)*eps_clr;
      phi[i][j][k]=phisurf+copysign(1.0,phi[i][j][k]-phisurf)*eps_clr;
    }
  }

  //clr.bnd_update();
  bdcond(clr);

#ifdef DEBUG
  std::cout<<"init:bnd_update \n";
#endif

#ifdef IB
  ib_ext_scalar(clr);
#ifdef DEBUG
  std::cout<<"init:ib_ext_scalar \n";
#endif
#endif

  clr.exchange_all();

  for_aijk(i,j,k){
    phi[i][j][k]=clr[i][j][k];
  }

  /* initial condition for fnode */
  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    scheme.f[i][j][k]=(phi[i-1][j-1][k-1]+phi[i  ][j-1][k-1]
                      +phi[i-1][j  ][k-1]+phi[i  ][j  ][k-1]
                      +phi[i-1][j-1][k  ]+phi[i  ][j-1][k  ]
                      +phi[i-1][j  ][k  ]+phi[i  ][j  ][k  ])/8.0;
  }}}

  scheme.bdcond_f(phi);

  for(int i=0; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    /* initial condition for line density */
    scheme.sigx[i][j][k]=(phi[i][j  ][k  ]+phi[i][j  ][k-1]
                         +phi[i][j-1][k  ]+phi[i][j-1][k-1])/4.0*dx;
  }}}

  scheme.bdcond_i(phi);

  for(int i=1; i<=ei()+1; i++){
  for(int j=0; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dy=phi.dyc(j);
    /* initial condition for line density */
    scheme.sigy[i][j][k]=(phi[i  ][j][k  ]+phi[i-1][j][k  ]
                         +phi[i  ][j][k-1]+phi[i-1][j][k-1])/4.0*dy;
  }}}

  scheme.bdcond_j(phi);

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=0; k<=ek()+1; k++){
    real dz=phi.dzc(k);
    /* initial condition for line density */
    scheme.sigz[i][j][k]=(phi[i  ][j  ][k]+phi[i-1][j  ][k]
                         +phi[i  ][j-1][k]+phi[i-1][j-1][k])/4.0*dz;
  }}}

  scheme.bdcond_k(phi);

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dy=phi.dyc(j);
    real dz=phi.dzc(k);
    /* initial condition for face density */
    sxyz[Comp::i()][i][j][k]=(phi[i][j][k]+phi[i-1][j][k])/2.0*dz*dy;
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dz=phi.dzc(k);
    /* initial condition for face density */
    sxyz[Comp::j()][i][j][k]=(phi[i][j][k]+phi[i][j-1][k])/2.0*dz*dx;
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dy=phi.dyc(j);
    /* initial condition for face density */
    sxyz[Comp::k()][i][j][k]=(phi[i][j][k]+phi[i][j][k-1])/2.0*dx*dy;
  }}}

  for_m(m)
    bdphiface(sxyz,m,phi);

#ifdef IB
  ib_bdcond(clr);
#endif
  scheme.f.exchange_all();
  scheme.sigx.exchange_all();
  scheme.sigy.exchange_all();
  scheme.sigz.exchange_all();
  sxyz.exchange_all();

#if 0
  color_minmax();
  if( maxval()>0.5 && minval()<0.5) {
    if(time->current_step()%nredist==0) {
      redist(false);
      boil::oout<<"# sharpening initial color function.\n";
    }
  }
#endif

#if 0
  plot_f("f.dat");
  plot_sigx("sigx.dat");
  plot_sigy("sigy.dat");
  plot_sigz("sigz.dat");
  plot_sxyz("sxyzi.dat",Comp::i());
  plot_sxyz("sxyzj.dat",Comp::j());
  plot_sxyz("sxyzk.dat",Comp::k());
#endif

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_init.cpp,v 1.10 2015/06/29 18:26:46 sato Exp $'/
+-----------------------------------------------------------------------------*/
