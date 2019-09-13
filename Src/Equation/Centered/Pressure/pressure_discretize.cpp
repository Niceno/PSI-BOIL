#include "pressure.h"

/***************************************************************************//**
*  Since discretization of pressure-Poisson equation is different from 
*  discretizations of other variables in a fact that diffusion-like
*  neighbouring coefficients depend on the inverse of physical property
*  (density), a specialized function is written here. 
*
*  The only thing useful from the parent is handling of boundary conditions.
*******************************************************************************/
void Pressure::discretize() {

  boil::timer.start("pressure discretize");

  Comp m;
  real a_w, a_e, a_s, a_n, a_b, a_t;
  real rhom, rhop;

  /*----------------------------+
  |  neighbouring coefficients  |
  +----------------------------*/

  /* coefficients in i direction (w and e) */
  m = Comp::u();
  for_ijk(i,j,k) {
    /* linear */
    rhom = fluid()->rho(m,i,  j,k);
    rhop = fluid()->rho(m,i+1,j,k);
    a_w = dSx(Sign::neg(),i,j,k);
    a_e = dSx(Sign::pos(),i,j,k);
    A.w[i][j][k] = a_w / dxw(i) / rhom;
    A.e[i][j][k] = a_e / dxe(i) / rhop;
  }

  /* coefficients in j direction (s and n) */
  m = Comp::v();
  for_ijk(i,j,k) {
    /* linear */
    rhom = fluid()->rho(m,i,j,  k);
    rhop = fluid()->rho(m,i,j+1,k);
    a_s = dSy(Sign::neg(),i,j,k);
    a_n = dSy(Sign::pos(),i,j,k);
    A.s[i][j][k] = a_s / dys(j) / rhom;
    A.n[i][j][k] = a_n / dyn(j) / rhop;
  }

  /* coefficients in k direction (b and t) */
  m = Comp::w();
  for_ijk(i,j,k) {
    /* linear */
    rhom = fluid()->rho(m,i,j,k);
    rhop = fluid()->rho(m,i,j,k+1);
    a_b = dSz(Sign::neg(),i,j,k);
    a_t = dSz(Sign::pos(),i,j,k);
    A.b[i][j][k] = a_b / dzb(k) / rhom;
    A.t[i][j][k] = a_t / dzt(k) / rhop;
  }

  Centered::create_system_bnd();

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) {
    for(int cc=0; cc<dom->ibody().nccells(); cc++) {
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

      /* both cells in fluid, but partially in solid */
      if( dom->ibody().cut_p(i,j,k) ) {
        if( dom->ibody().cut_p(i-1,j,k) ) A.w[i][j][k] *= dom->ibody().fSw(i,j,k);
        if( dom->ibody().cut_p(i+1,j,k) ) A.e[i][j][k] *= dom->ibody().fSe(i,j,k);
        if( dom->ibody().cut_p(i,j-1,k) ) A.s[i][j][k] *= dom->ibody().fSs(i,j,k);
        if( dom->ibody().cut_p(i,j+1,k) ) A.n[i][j][k] *= dom->ibody().fSn(i,j,k);
        if( dom->ibody().cut_p(i,j,k-1) ) A.b[i][j][k] *= dom->ibody().fSb(i,j,k);
        if( dom->ibody().cut_p(i,j,k+1) ) A.t[i][j][k] *= dom->ibody().fSt(i,j,k);
      }

      /* some cells are inside the immersed body 
         (there seems to be  too many checks, each
          face is checked twice, but it is needed) */
      if( dom->ibody().off_p(i,j,k) ) {
        A.w[i][j][k] = 0.0; A.e[i-1][j][k] = 0.0;
        A.e[i][j][k] = 0.0; A.w[i+1][j][k] = 0.0;
        A.s[i][j][k] = 0.0; A.n[i][j-1][k] = 0.0; 
        A.n[i][j][k] = 0.0; A.s[i][j+1][k] = 0.0;
        A.b[i][j][k] = 0.0; A.t[i][j][k-1] = 0.0; 
        A.t[i][j][k] = 0.0; A.b[i][j][k+1] = 0.0;
      }

      if(dom->ibody().off_p(i-1,j,k)) {A.w[i][j][k] = 0.0; A.e[i-1][j][k] = 0.0;}
      if(dom->ibody().off_p(i+1,j,k)) {A.e[i][j][k] = 0.0; A.w[i+1][j][k] = 0.0;}
      if(dom->ibody().off_p(i,j-1,k)) {A.s[i][j][k] = 0.0; A.n[i][j-1][k] = 0.0;}
      if(dom->ibody().off_p(i,j+1,k)) {A.n[i][j][k] = 0.0; A.s[i][j+1][k] = 0.0;}
      if(dom->ibody().off_p(i,j,k-1)) {A.b[i][j][k] = 0.0; A.t[i][j][k-1] = 0.0;}
      if(dom->ibody().off_p(i,j,k+1)) {A.t[i][j][k] = 0.0; A.b[i][j][k+1] = 0.0;}
    }

  } /* is there an immersed body */

  /*----------------------+
  |  central coefficient  |
  +----------------------*/
  for_ijk(i,j,k) A.c[i][j][k] = A.w[i][j][k] + A.e[i][j][k] 
                              + A.s[i][j][k] + A.n[i][j][k]
                              + A.b[i][j][k] + A.t[i][j][k];

  for_ijk(i,j,k) A.ci[i][j][k] = 1.0 / A.c[i][j][k];

  /*-------------------------------+
  |  a "touch" from immersed body  | -> this should be necessary only for A.c 
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) 
    for_ijk(i,j,k) 
      if( dom->ibody().off_p(i,j,k) ) {
        A.c[i][j][k]  = boil::pico;
        A.w[i][j][k]  = 0.0;
        A.e[i][j][k]  = 0.0;
        A.s[i][j][k]  = 0.0;
        A.n[i][j][k]  = 0.0;
        A.b[i][j][k]  = 0.0;
        A.t[i][j][k]  = 0.0;
        A.ci[i][j][k] = 1.0;
    }

  boil::timer.stop("pressure discretize");
}
