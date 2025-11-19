#include "vofaxisym.h"

//#define CRAY

/***************************************************************************//**
 \brief Calculate the line constant alpha for given configuration. Moreover,
        some threshold value functions are included.
        small = mmx<mmz, negative = nx<0 and vice versa
        source: derived by me (= Lubomir)
*******************************************************************************/
real VOFaxisym::phi_max_small_pos(const real mmx,const real mmz,const real eta0) {
  real xi = xi_small_pos_trapezoid(0.5,mmx,mmz);
  real c = calc_v(0.5,mmx,0.0,mmz);
  real K = (eta0+xi)/(eta0+0.5);
  return c*K;
}

real VOFaxisym::phi_max_small_neg(const real mmx,const real mmz,const real eta0) {
  real xi = xi_small_neg_trapezoid(0.5,mmx,mmz);
  real c = calc_v(0.5,mmx,0.0,mmz);
  real K = (eta0+xi)/(eta0+0.5);
  return c*K;
}

real VOFaxisym::phi_max_large_pos(const real mmx,const real mmz,const real eta0) {
  real xi = xi_large_pos_trapezoid(0.5,mmx,mmz);
  real c = calc_v(0.5,mmx,0.0,mmz);
  real K = (eta0+xi)/(eta0+0.5);
  return c*K;
}

real VOFaxisym::phi_max_large_neg(const real mmx,const real mmz,const real eta0) {
  real xi = xi_large_neg_trapezoid(0.5,mmx,mmz);
  real c = calc_v(0.5,mmx,0.0,mmz);
  real K = (eta0+xi)/(eta0+0.5);
  return c*K;
}

real VOFaxisym::phi_tr_small_pos(const real mmx,const real mmz,const real eta0) {
  return mmx/2./mmz*(eta0+1./3.)/(eta0+0.5);
}

real VOFaxisym::phi_tr_small_neg(const real mmx,const real mmz,const real eta0) {
  return mmx/2./mmz*(eta0+2./3.)/(eta0+0.5);
}

real VOFaxisym::phi_tr_large_pos(const real mmx,const real mmz,const real eta0) {
#ifndef CRAY
   return mmz/2./mmx*(eta0+mmz/3./mmx)/(eta0+0.5);
#else
   return mmx>0.0 ? mmz/2./mmx*(eta0+mmz/3./mmx)/(eta0+0.5) : 0.0;
#endif
}

real VOFaxisym::phi_tr_large_neg(const real mmx,const real mmz,const real eta0) {
#ifndef CRAY
   return mmz/2./mmx*(eta0+1.-mmz/3./mmx)/(eta0+0.5);
#else
   return mmx>0.0 ? mmz/2./mmx*(eta0+1.-mmz/3./mmx)/(eta0+0.5) : 0.0;
#endif
}

real VOFaxisym::phi_crit_pos(const real mmx,const real mmz,const real eta0) {
  return 2.*mmx/3./mmz * pow(eta0,3.)/(eta0+0.5);
}

real VOFaxisym::phi_crit_neg(const real mmx,const real mmz,const real eta0) {
  return 2.*mmx/3./mmz * pow(eta0+1.,3.)/(eta0+0.5);
}

real VOFaxisym::alp_pos_triangle_subcrit(const real vf,const real mmx,const real mmz,
                                         const real eta0,const real phi_crit) {
  real M = 3.*mmx*mmx*mmz*(eta0+0.5);
  real angle = asin(1.-2.*vf/phi_crit);
  return pow(M*phi_crit/2.,1./3.) * (-1. + 2.*cos(angle/3.+boil::pi/6.));
}

real VOFaxisym::alp_pos_triangle_supercrit(const real vf,const real mmx,const real mmz,
                                           const real eta0,const real phi_crit) {
  real M = 3.*mmx*mmx*mmz*(eta0+0.5);
  real coef = 1./3.;
  return pow(M/2.,coef) * (-pow(phi_crit,coef)
                           +pow(2.*vf-phi_crit+2.*sqrt(vf*(vf-phi_crit)),coef)
                           +pow(2.*vf-phi_crit-2.*sqrt(vf*(vf-phi_crit)),coef) );
}

real VOFaxisym::alp_neg_triangle(const real vf,const real mmx,const real mmz,const real eta0) {
  real M = 3.*mmx*mmx*mmz*(eta0+0.5);
  real phi_crit = phi_crit_neg(mmx,mmz,eta0);
  real angle = asin(1.-2.*vf/phi_crit);
  return -pow(M*phi_crit/2.,1./3.) * (-1.+2.*sin(angle/3.));
}

real VOFaxisym::alp_small_pos_trapezoid(const real vf,const real mmx,
                                        const real mmz,const real eta0) {
  return mmz*vf + mmx/2.*(eta0+2./3.)/(eta0+0.5);
}

real VOFaxisym::alp_small_neg_trapezoid(const real vf,const real mmx,
                                        const real mmz,const real eta0) {
  return mmz*vf + mmx/2.*(eta0+1./3.)/(eta0+0.5);
}

real VOFaxisym::alp_large_pos_trapezoid(const real vf,const real mmx,
                                        const real mmz,const real eta0) {
  return (mmz-2.*mmx*eta0+sqrt(4.*mmx*mmx*eta0*eta0
         +8.*mmx*mmx*(eta0+0.5)*vf-mmz*mmz/3.      ))/2.;
}

real VOFaxisym::alp_large_neg_trapezoid(const real vf,const real mmx,
                                        const real mmz,const real eta0) {
  return (mmz+2.*mmx*(eta0+1.)-sqrt(4.*mmx*mmx*pow(eta0+1.,2.)
         -8.*mmx*mmx*(eta0+0.5)*vf-mmz*mmz/3.                 ))/2.;
}
