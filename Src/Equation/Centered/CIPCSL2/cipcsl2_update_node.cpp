#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::update_node(Scalar & sca) {
/***************************************************************************//**
*  \brief update node, edge, face and cell-centered value.
*         input : sca
*         original color function: clr
*         tempolary: fn
*******************************************************************************/

  /*--------------------------+
  |  update scheme variables  |
  +--------------------------*/

#ifdef IB
  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  ib_ext_scalar(sca);
  ib_bdcond(sca);
#endif


  /* calculate delt clr */
  for_aijk(i,j,k)
    fn[i][j][k]=sca[i][j][k]-clr[i][j][k];

  /* update fnode */
  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real tmp = scheme.f[i][j][k] + (fn[i-1][j-1][k-1]+fn[i  ][j-1][k-1]
                                   +fn[i-1][j  ][k-1]+fn[i  ][j  ][k-1]
                                   +fn[i-1][j-1][k  ]+fn[i  ][j-1][k  ]
                                   +fn[i-1][j  ][k  ]+fn[i  ][j  ][k  ])/8.0;
    scheme.f[i][j][k]=max(0.0,min(1.0,tmp));
  }}}
  scheme.bdcond_f(sca);
  scheme.f.exchange_all();

  for(int i=0; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real tmp = scheme.sigx[i][j][k] + (fn[i][j  ][k  ]+fn[i][j  ][k-1]
                                      +fn[i][j-1][k  ]+fn[i][j-1][k-1])/4.0*dx;
    scheme.sigx[i][j][k] = max(0.0,min(dx,tmp));
  }}}
  scheme.bdcond_i(sca);
  scheme.sigx.exchange_all();


  for(int i=1; i<=ei()+1; i++){
  for(int j=0; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dy=phi.dyc(j);
    real tmp = scheme.sigy[i][j][k] + (fn[i  ][j][k  ]+fn[i-1][j][k  ]
                                      +fn[i  ][j][k-1]+fn[i-1][j][k-1])/4.0*dy;
    scheme.sigy[i][j][k] = max(0.0,min(dy,tmp)); 
  }}}
  scheme.bdcond_j(sca);
  scheme.sigy.exchange_all();


  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=0; k<=ek()+1; k++){
    real dz=phi.dzc(k);
    real tmp = scheme.sigz[i][j][k] + (fn[i  ][j  ][k]+fn[i-1][j  ][k]
                                      +fn[i  ][j-1][k]+fn[i-1][j-1][k])/4.0*dz;
    scheme.sigz[i][j][k] = max(0.0,min(dz,tmp));
  }}}
  scheme.bdcond_k(sca);
  scheme.sigz.exchange_all();

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dy=phi.dyc(j);
    real dz=phi.dzc(k);
    real tmp = sxyz[Comp::i()][i][j][k] + (fn[i][j][k]+fn[i-1][j][k])/2.0*dz*dy;
    sxyz[Comp::i()][i][j][k] = max(0.0,min(dz*dy,tmp));
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dz=phi.dzc(k);
    real tmp = sxyz[Comp::j()][i][j][k] + (fn[i][j][k]+fn[i][j-1][k])/2.0*dz*dx;
    sxyz[Comp::j()][i][j][k] = max(0.0,min(dz*dx,tmp));
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dy=phi.dyc(j);
    real tmp = sxyz[Comp::k()][i][j][k] + (fn[i][j][k]+fn[i][j][k-1])/2.0*dx*dy;
    sxyz[Comp::k()][i][j][k] = max(0.0,min(dx*dy,tmp));
  }}}
  for_m(m)
    bdphiface(sxyz,m,sca);
  sxyz.exchange_all();

  /*-------------+
  |  update clr  |
  +-------------*/
  /* calculate delt clr */
  for_aijk(i,j,k)
    clr[i][j][k]=sca[i][j][k];

#ifdef IB
  /*--------------------+
  |  immersed boundary  |
  +--------------------*/
  ib_ext_scalar(clr);
  ib_bdcond(clr); 
#endif

  /* ancillary functions */
  ancillary();

  return;
}
