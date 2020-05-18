#include "cavitypressure.h"

/***************************************************************************//**
*  Discretization near the phasic interface.
*  Inherits structure from EnthalpyFD.
*******************************************************************************/
void CavityPressure::diff_matrix(const int i, const int j, const int k) {

  /* reset */
  A.c[i][j][k] = 0.0;
#if 0

  /* auxilliary parameters */
  const real mult = phi.dV(i,j,k)/fluid()->rho(i,j,k);
  real xm,xp,aflagm,aflagp;
  aflagm=aflagp=1.0;

  /* coefficients in i direction (w and e) */
  if(!interface(Sign::neg(),Comp::i(),i,j,k)) {
    xm = phi.dxw(i);
  } else {
    real ps;
    xm = std::max(epsl*phi.dxw(i),distance_int_x(Sign::neg(),i,j,k,ps));
    aflagm=0.0;
  }
  if(!interface(Sign::pos(),Comp::i(),i,j,k)){
    xp=phi.dxe(i);
  } else {
    real ps;
    xp = std::max(epsl*phi.dxe(i),distance_int_x(Sign::pos(),i,j,k,ps));
    aflagp=0.0;
  }
  real cxm = coef_x_m(xm,xp,phi.xc(i));
  real cxp = coef_x_p(xm,xp,phi.xc(i));

  A.w[i][j][k] =  mult * cxm * aflagm;
  A.c[i][j][k] += mult * (cxm+cxp);
  A.e[i][j][k] =  mult * cxp * aflagp;

  /* coefficients in j direction (s and n) */
  real ym,yp;
  aflagm=aflagp=1.0;

  if(!interface(Sign::neg(),Comp::j(),i,j,k)){
    ym=phi.dys(j);
  } else {
    real ps;
    ym = std::max(epsl*phi.dys(j),distance_int_y(Sign::neg(),i,j,k,ps));
    aflagm=0.0;
  }
  if(!interface(Sign::pos(),Comp::j(),i,j,k)){
    yp=phi.dyn(j);
  } else {
    real ps;
    yp = std::max(epsl*phi.dyn(j),distance_int_y(Sign::pos(),i,j,k,ps));
    aflagp=0.0;
  }
  real cym = coef_y_m(ym,yp,phi.yc(j));
  real cyp = coef_y_p(ym,yp,phi.yc(j));

  A.s[i][j][k] =  mult * cym * aflagm;
  A.c[i][j][k] += mult * (cym+cyp);
  A.n[i][j][k] =  mult * cyp * aflagp;

  /* coefficients in k direction (b and t) */
  real zm,zp;
  aflagm=aflagp=1.0;

  if(!interface(Sign::neg(),Comp::k(),i,j,k)){
    zm=phi.dzb(k);
  } else {
    real ps;
    zm = std::max(epsl*phi.dzb(k),distance_int_z(Sign::neg(),i,j,k,ps));
    aflagm=0.0;
  }
  if(!interface(Sign::pos(),Comp::k(),i,j,k)){
    zp=phi.dzt(k);
  } else {
    real ps;
    zp = std::max(epsl*phi.dzt(k),distance_int_z(Sign::pos(),i,j,k,ps));
    aflagp=0.0;
  }
  real czm = coef_z_m(zm,zp,phi.zc(k));
  real czp = coef_z_p(zm,zp,phi.zc(k));

  A.b[i][j][k] =  mult * czm * aflagm;
  A.c[i][j][k] += mult * (czm+czp);
  A.t[i][j][k] =  mult * czp * aflagp;

#endif
  return;
}
