#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
real PhaseChangeVOF::gradtx8(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in x-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  Comp mcomp = Comp::i();
  int m1, m2;
  real dxm1, dxm2;
  real tm0 = tpr[i][j][k];
  real tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    m2 = 2;
    dxm1 = clr.dxw(i+m1);
    tm1  = tpr[i+m1][j][k];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (i==ei())
      if (!tpr.bc().type(Dir::imax(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::imax())) {
        return (tm1-tm0)/dxm1;
      }

    real temp = bndtpr[mcomp][i+m2][j][k];
    if (!boil::realistic(temp)) {
      dxm2 = clr.dxc(i)/2.0;
      tm2  = bndtpr[mcomp][i+m1][j][k];
    } else {
      dxm2 = clr.dxw(i+m1) + clr.dxc(i+m1)/2.0;
      tm2  = temp;
    }
  } else {
    m1 = -1;
    m2 = -2;
    dxm1 = -clr.dxe(i+m1);
    tm1  = tpr[i+m1][j][k];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (i==si())
      if (!tpr.bc().type(Dir::imin(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::imin())) {
        return (tm1-tm0)/dxm1;
      }

    real temp = bndtpr[mcomp][i+m1][j][k];
    if (!boil::realistic(temp)) {
      dxm2 = -clr.dxc(i)/2.0;
      tm2  = bndtpr[mcomp][i   ][j][k];
    } else {
      dxm2 = -clr.dxe(i+m1) - clr.dxc(i+m1)/2.0;
      tm2  = temp;
    }
  }

  return grad_2nd(tm0,tm1,tm2,dxm1,dxm2);
}

/******************************************************************************/
real PhaseChangeVOF::gradty8(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in y-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  Comp mcomp = Comp::j();
  int m1, m2;
  real dym1, dym2;
  real tm0 = tpr[i][j][k];
  real tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    m2 = 2;
    dym1 = clr.dys(j+m1);
    tm1  = tpr[i][j+m1][k];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (j==ej())
      if (!tpr.bc().type(Dir::jmax(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::jmax())) {
        return (tm1-tm0)/dym1;
      }
    
    real temp = bndtpr[mcomp][i][j+m2][k];
    if (!boil::realistic(temp)) {
      dym2 = clr.dyc(j)/2.0;
      tm2  = bndtpr[mcomp][i][j+m1][k];
    } else {
      dym2 = clr.dys(j+m1) + clr.dyc(j+m1)/2.0;
      tm2  = temp;
    }
  } else {
    m1 = -1;
    m2 = -2;
    dym1 = -clr.dyn(j+m1);
    tm1  = tpr[i][j+m1][k];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (j==sj())
      if (!tpr.bc().type(Dir::jmin(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::jmin())) {
        return (tm1-tm0)/dym1;
      }

    real temp = bndtpr[mcomp][i][j+m1][k];
    if (!boil::realistic(temp)) {
      dym2 = -clr.dyc(j)/2.0;
      tm2  = bndtpr[mcomp][i][j   ][k];
    } else {
      dym2 = -clr.dyn(j+m1) - clr.dyc(j+m1)/2.0;
      tm2  = temp;
    }
  }

  return grad_2nd(tm0,tm1,tm2,dym1,dym2);
}

/******************************************************************************/
real PhaseChangeVOF::gradtz8(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in z-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  Comp mcomp = Comp::k();
  int m1, m2;
  real dzm1, dzm2;
  real tm0 = tpr[i][j][k];
  real tm1, tm2;
  if (dir < 0) {
    m1 = 1;
    m2 = 2;
    dzm1 = clr.dzb(k+m1);
    tm1  = tpr[i][j][k+m1];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (k==ek())
      if (!tpr.bc().type(Dir::kmax(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::kmax())) {
        return (tm1-tm0)/dzm1;
      }

    real temp = bndtpr[mcomp][i][j][k+m2];
    if (!boil::realistic(temp)) {
      dzm2 = clr.dzc(k)/2.0;
      tm2  = bndtpr[mcomp][i][j][k+m1];
    } else {
      dzm2 = clr.dzb(k+m1) + clr.dzc(k+m1)/2.0;
      tm2  = temp;
    }
  } else {
    m1 = -1;
    m2 = -2;
    dzm1 = -clr.dzt(k+m1);
    tm1  = tpr[i][j][k+m1];

    /*at non-periodic boundaries, not enough points exist for 2nd order scheme*/
    if (k==sk())
      if (!tpr.bc().type(Dir::kmin(),BndType::periodic())&&!tpr.bc().type_decomp(Dir::kmin())) {
        return (tm1-tm0)/dzm1;
      }

    real temp = bndtpr[mcomp][i][j][k+m1];
    if (!boil::realistic(temp)) {
      dzm2 = -clr.dzc(k)/2.0;
      tm2  = bndtpr[mcomp][i][j][k   ];
    } else {
      dzm2 = -clr.dzt(k+m1) - clr.dzc(k+m1)/2.0;
      tm2  = temp;
    }
  }

  return grad_2nd(tm0,tm1,tm2,dzm1,dzm2);
}
