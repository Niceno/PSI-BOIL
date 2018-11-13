#include "phasechange.h"
using namespace std;

/******************************************************************************/
void PhaseChange::prepare_gradt8() {
/***************************************************************************//*** 
*  \brief prepare gradt in the domain
*  bndtpr is positive for liquid and negative for vapor
*******************************************************************************/
  Comp m;
 
  m = Comp::i();
  for(int i=si(); i<=ei()+1; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = tpr[i][j][k]*tpr.dxc(i-1) + tpr[i-1][j][k]*tpr.dxc(i);
    temp /= (tpr.dxc(i) + tpr.dxc(i-1));

    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::j();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()+1; j++)
  for(int k=sk(); k<=ek()  ; k++) {
    real temp = tpr[i][j][k]*tpr.dyc(j-1) + tpr[i][j-1][k]*tpr.dyc(j);
    temp /= (tpr.dyc(j) + tpr.dyc(j-1));

    bndtpr[m][i][j][k] = temp;
  }

  m = Comp::k();
  for(int i=si(); i<=ei()  ; i++)
  for(int j=sj(); j<=ej()  ; j++)
  for(int k=sk(); k<=ek()+1; k++) {
    real temp = tpr[i][j][k]*tpr.dzc(k-1) + tpr[i][j][k-1]*tpr.dzc(k);
    temp /= (tpr.dzc(k) + tpr.dzc(k-1));

    bndtpr[m][i][j][k] = temp;
  }

  bndtpr.exchange_all();
}

/******************************************************************************/
real PhaseChange::gradtx8(const int dir, const int i, const int j, const int k){
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i+m2,j,k)-phisurf)<=0.0) {
      dxm2 = clr.dxc(i)/2.0;
      tm2  = bndtpr[mcomp][i+m1][j][k];
    } else {
      dxm2 = clr.dxw(i+m1) + clr.dxc(i+m1)/2.0;
      tm2  = bndtpr[mcomp][i+m2][j][k];
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i+m1,j,k)-phisurf)<=0.0) {
      dxm2 = -clr.dxc(i)/2.0;
      tm2  = bndtpr[mcomp][i   ][j][k];
    } else {
      dxm2 = -clr.dxe(i+m1) - clr.dxc(i+m1)/2.0;
      tm2  = bndtpr[mcomp][i+m1][j][k];
    }
  }

  return gradt8(tm0,tm1,tm2,dxm1,dxm2);
}

/******************************************************************************/
real PhaseChange::gradty8(const int dir, const int i, const int j, const int k){
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i,j+m2,k)-phisurf)<=0.0) {
      dym2 = clr.dyc(j)/2.0;
      tm2  = bndtpr[mcomp][i][j+m1][k];
    } else {
      dym2 = clr.dys(j+m1) + clr.dyc(j+m1)/2.0;
      tm2  = bndtpr[mcomp][i][j+m2][k];
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i,j+m1,k)-phisurf)<=0.0) {
      dym2 = -clr.dyc(j)/2.0;
      tm2  = bndtpr[mcomp][i][j   ][k];
    } else {
      dym2 = -clr.dyn(j+m1) - clr.dyc(j+m1)/2.0;
      tm2  = bndtpr[mcomp][i][j+m1][k];
    }
  }

  return gradt8(tm0,tm1,tm2,dym1,dym2);
}

/******************************************************************************/
real PhaseChange::gradtz8(const int dir, const int i, const int j, const int k){
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i,j,k+m2)-phisurf)<=0.0) {
      dzm2 = clr.dzc(k)/2.0;
      tm2  = bndtpr[mcomp][i][j][k+m1];
    } else {
      dzm2 = clr.dzb(k+m1) + clr.dzc(k+m1)/2.0;
      tm2  = bndtpr[mcomp][i][j][k+m2];
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

    if ((clr[i][j][k]-phisurf) *
        (surface_color(mcomp,i,j,k+m1)-phisurf)<=0.0) {
      dzm2 = -clr.dzc(k)/2.0;
      tm2  = bndtpr[mcomp][i][j][k   ];
    } else {
      dzm2 = -clr.dzt(k+m1) - clr.dzc(k+m1)/2.0;
      tm2  = bndtpr[mcomp][i][j][k+m1];
    }
  }

  return gradt8(tm0,tm1,tm2,dzm1,dzm2);
}

/******************************************************************************/
real PhaseChange::gradt8(const real tm0, const real tm1, const real tm2,
			 const real dm1, const real dm2){
/***************************************************************************//**
*  \brief calculate second order accurate grad(tpr) 
*******************************************************************************/
  real denom = dm1*dm2*(dm2-dm1);

  return dm2*dm2*(tm1-tm0)/denom - dm1*dm1*(tm2-tm0)/denom;
}

