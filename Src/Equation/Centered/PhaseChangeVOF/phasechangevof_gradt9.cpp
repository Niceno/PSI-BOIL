#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
real PhaseChangeVOF::gradtx9(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in x-direction
*         dir= 1: grad(tpr) second order central
*         dir=-1: grad(tpr)  second order cenral
*******************************************************************************/
  int m1;
  real dxm1, dxm2;
  real tm0, tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    dxm1 = distance_x(i,j,k,dir,tm0);
    tm1  = tpr[i][j][k];

    /* is the interface too close to the cell centre? */
    if(dxm1/clr.dxc(i) < epsl)
      return gradtx8(dir,i,j,k);
  
    dxm2 = dxm1 + clr.dxe(i); 
    tm2  = tpr[i+m1][j][k];   

    /* cell-centred gradient */
    return grad_2nd(tm1,tm0,tm2,-dxm1,dxm2-dxm1);
  } else {
    m1 = -1;
    dxm1 = -distance_x(i,j,k,dir,tm0);
    tm1  = tpr[i][j][k];

    /* is the interface too close to the cell centre? */
    if(fabs(dxm1/clr.dxc(i)) < epsl)
      return gradtx8(dir,i,j,k);

    dxm2 = dxm1 - clr.dxw(i);
    tm2  = tpr[i+m1][j][k];   

    /* cell-centred gradient */
    return grad_2nd(tm1,tm2,tm0,dxm2-dxm1,-dxm1);
  }

  /* upwind/downwind at interface */
  //return grad_2nd(tm0,tm1,tm2,dxm1,dxm2);
}

/******************************************************************************/
real PhaseChangeVOF::gradty9(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in y-direction
*         dir= 1: grad(tpr) second order central
*         dir=-1: grad(tpr)  second order central
*******************************************************************************/
  int m1;
  real dym1, dym2;
  real tm0, tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    dym1 = distance_y(i,j,k,dir,tm0);

    /* is the interface too close to the cell centre? */
    if(dym1/clr.dyc(j) < epsl)
      return gradty8(dir,i,j,k);

    tm1  = tpr[i][j][k];

    dym2 = dym1 + clr.dyn(j);   
    tm2  = tpr[i][j+m1][k];

    /* cell-centred gradient */
    return grad_2nd(tm1,tm0,tm2,-dym1,dym2-dym1);
  } else {
    m1 = -1;
    dym1 = -distance_y(i,j,k,dir,tm0);
    tm1  = tpr[i][j][k];

    /* is the interface too close to the cell centre? */
    if(fabs(dym1/clr.dyc(j)) < epsl)
      return gradty8(dir,i,j,k);

    dym2 = dym1 - clr.dys(j);
    tm2  = tpr[i][j+m1][k];

    /* cell-centred gradient */
    return grad_2nd(tm1,tm2,tm0,dym2-dym1,-dym1);
  }

  /* upwind/downwind at interface */
  //return grad_2nd(tm0,tm1,tm2,dym1,dym2);
}

/******************************************************************************/
real PhaseChangeVOF::gradtz9(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in z-direction
*         dir= 1: grad(tpr) second order central
*         dir=-1: grad(tpr)  second order central
*******************************************************************************/
  int m1;
  real dzm1, dzm2;
  real tm0, tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    dzm1 = distance_z(i,j,k,dir,tm0);
    tm1  = tpr[i][j][k];

    /* is the interface too close to the cell centre? */
    if(dzm1/clr.dzc(k) < epsl)
      return gradtz8(dir,i,j,k);

    dzm2 = dzm1 + clr.dzt(k);
    tm2  = tpr[i][j][k+m1];

    /* cell-centred gradient */
    return grad_2nd(tm1,tm0,tm2,-dzm1,dzm2-dzm1);
  } else {
    m1 = -1;
    dzm1 = -distance_z(i,j,k,dir,tm0);
    tm1  = tpr[i][j][k];

    /* is the interface too close to the cell centre? */
    if(fabs(dzm1/clr.dzc(k)) < epsl)
      return gradtz8(dir,i,j,k);

    dzm2 = dzm1 - clr.dzb(k);
    tm2  = tpr[i][j][k+m1];

    /* cell-centred gradient */
    return grad_2nd(tm1,tm2,tm0,dzm2-dzm1,-dzm1);
  }

  /* upwind/downwind at interface */
  //return grad_2nd(tm0,tm1,tm2,dzm1,dzm2);
}
