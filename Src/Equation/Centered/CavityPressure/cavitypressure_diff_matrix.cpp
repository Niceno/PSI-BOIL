#include "cavitypressure.h"

/***************************************************************************//**
*  Discretization of system matrix.
*  Inherits structure from EnthalpyFD.
*******************************************************************************/
void CavityPressure::diff_matrix(const int i, const int j, const int k) {

  /* reset */
  A.c[i][j][k]  = 0.0;

  /* auxilliary parameters */
  const real mult = phi.dV(i,j,k)/fluid()->rho(i,j,k);
  real xm,xp;
  real pm,pp;
  int aflagm,aflagp;
  aflagm=aflagp=1;

  /* coefficients in i direction (w and e) */
  if(!interface(Sign::neg(),Comp::i(),i,j,k)) {
    xm = phi.dxw(i);
    pm = phi[i-1][j][k];
  } else {
    xm = distance_int_x(Sign::neg(),i,j,k,pm);
    aflagm=0;
  }
  if(!interface(Sign::pos(),Comp::i(),i,j,k)){
    xp = phi.dxe(i);
    pp = phi[i+1][j][k];
  } else {
    xp = distance_int_x(Sign::pos(),i,j,k,pp);
    aflagp=0;
  }
  real cxm = coef_x_m(xm,xp,phi.xc(i));
  real cxp = coef_x_p(xm,xp,phi.xc(i));

  A.w[i][j][k] =  mult * cxm * aflagm;
  A.c[i][j][k] += mult * (cxm+cxp);
  A.e[i][j][k] =  mult * cxp * aflagp;

  /* explicit term */
  fold[i][j][k] += mult*cxm*(1-aflagm)*pm
                 + mult*cxp*(1-aflagp)*pp;

  /* coefficients in j direction (s and n) */
  real ym,yp;
  aflagm=aflagp=1;

  if(!interface(Sign::neg(),Comp::j(),i,j,k)){
    ym = phi.dys(j);
    pm = phi[i][j-1][k];
  } else {
    ym = distance_int_y(Sign::neg(),i,j,k,pm);
    aflagm=0;
  }
  if(!interface(Sign::pos(),Comp::j(),i,j,k)){
    yp = phi.dyn(j);
    pp = phi[i][j+1][k];
  } else {
    yp = distance_int_y(Sign::pos(),i,j,k,pp);
    aflagp=0;
  }
  real cym = coef_y_m(ym,yp,phi.yc(j));
  real cyp = coef_y_p(ym,yp,phi.yc(j));

  A.s[i][j][k] =  mult * cym * aflagm;
  A.c[i][j][k] += mult * (cym+cyp);
  A.n[i][j][k] =  mult * cyp * aflagp;

  /* explicit term */
  fold[i][j][k] += mult*cym*(1-aflagm)*pm
                 + mult*cyp*(1-aflagp)*pp;

  /* coefficients in k direction (b and t) */
  real zm,zp;
  aflagm=aflagp=1;

  if(!interface(Sign::neg(),Comp::k(),i,j,k)){
    zm = phi.dzb(k);
    pm = phi[i][j][k-1];
  } else {
    zm = distance_int_z(Sign::neg(),i,j,k,pm);
    aflagm=0;
  }
  if(!interface(Sign::pos(),Comp::k(),i,j,k)){
    zp = phi.dzt(k);
    pp = phi[i][j][k+1];
  } else {
    zp = distance_int_z(Sign::pos(),i,j,k,pp);
    aflagp=0;
  }
  real czm = coef_z_m(zm,zp,phi.zc(k));
  real czp = coef_z_p(zm,zp,phi.zc(k));

  A.b[i][j][k] =  mult * czm * aflagm;
  A.c[i][j][k] += mult * (czm+czp);
  A.t[i][j][k] =  mult * czp * aflagp;

  /* explicit term */
  fold[i][j][k] += mult*czm*(1-aflagm)*pm
                 + mult*czp*(1-aflagp)*pp;

  return;
}
