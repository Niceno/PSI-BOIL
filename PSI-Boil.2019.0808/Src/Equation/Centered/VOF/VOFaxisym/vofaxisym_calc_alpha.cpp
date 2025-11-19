#include "vofaxisym.h"

#if 0 /* discarded due to the way elvira works */
/******************************************************************************/
real VOFaxisym::calc_alpha(const real v,
                           const real vma, const real vmb, const real vmc)  const {
/***************************************************************************//**
 \brief Solve the point-inverse 2D Cartesian problem, i.e. calculate alp(phi,n),
        this is a virtual function wrt 3D VOF calc_v; these two should basically
        produce the same result (2D = special case 3D for Cartesian);
        source: paper by Scardovelli and Zaleski(2000)
*******************************************************************************/

  real w, v1, vm1, vm3; 
  /* vm2 = vmb = 0! This is 2D xz! */
  real alpha;

  w     = boil::minr(v, 1.0-v);  
  vm1   = boil::minr(vma,vmc);
  vm3   = 1.0-vm1;
  v1    = 0.5*vm1/vm3;

  if(w<v1) {
    alpha = sqrt(2.*vm1*vm3*w);
  } else {
    alpha = vm3*(w+v1);
  }

  if(v > 0.5){
      alpha = 1.0 - alpha;
  }
  
  return alpha;
}
#endif

/******************************************************************************/
real VOFaxisym::calc_alpha_axisymmetric(const real nnnx, const real v,
                                        const real eta0) {
/***************************************************************************//**
 \brief Solve the point-inverse axisymmetric problem, i.e. calc alp(phi,n,eta),
        source: derived by me (= Lubomir)
*******************************************************************************/
  if(v<boil::pico)
    return 0.0;
  if(v-1.>-boil::pico)
    return 1.0;

  real nnx = nnnx;
  real vf = v;

  /* nnx is already normalized by L1 norm */
  real mmx = fabs(nnx);
  real mmz = 1.0-mmx;

  /* complementary problem */
  real vf_max(0.5);

  if(mmx<=mmz) {
    if(nnx>=0.) {
      vf_max = phi_max_small_pos(mmx,mmz,eta0);
    } else {
      vf_max = phi_max_small_neg(mmx,mmz,eta0);
    }
  } else {
    if(nnx>=0.) {
      vf_max = phi_max_large_pos(mmx,mmz,eta0);
    } else {
      vf_max = phi_max_large_neg(mmx,mmz,eta0);
    }
  }

  if(v>vf_max) {
    vf = 1. - vf;
    nnx = -nnx;
  }

  real alp = 1.0;

  /* axisymmetric solution */
  if(mmx <= mmz) {
    if(nnx >= 0.) {
      if(vf < phi_tr_small_pos(mmx,mmz,eta0)) {
        real phi_crit = phi_crit_pos(mmx,mmz,eta0);
        if(vf < phi_crit) {
          alp = alp_pos_triangle_subcrit(vf,mmx,mmz,eta0,phi_crit);
        } else {
          alp = alp_pos_triangle_supercrit(vf,mmx,mmz,eta0,phi_crit);
        }
      } else {
        alp = alp_small_pos_trapezoid(vf,mmx,mmz,eta0);
      }
    } else {
      if(vf < phi_tr_small_neg(mmx,mmz,eta0)) {
        alp = alp_neg_triangle(vf,mmx,mmz,eta0);
      } else {
        alp = alp_small_neg_trapezoid(vf,mmx,mmz,eta0);
      }
    }
  } else {
    if(nnx >= 0) {
      if(vf < phi_tr_large_pos(mmx,mmz,eta0)) {
        real phi_crit = phi_crit_pos(mmx,mmz,eta0);
        if(vf < phi_crit) {
          alp = alp_pos_triangle_subcrit(vf,mmx,mmz,eta0,phi_crit);
        } else {
          alp = alp_pos_triangle_supercrit(vf,mmx,mmz,eta0,phi_crit);
        }
      } else {
        alp = alp_large_pos_trapezoid(vf,mmx,mmz,eta0);
      }
    } else {
      if(vf < phi_tr_large_neg(mmx,mmz,eta0)) {
        alp = alp_neg_triangle(vf,mmx,mmz,eta0);
      } else {
        alp = alp_large_neg_trapezoid(vf,mmx,mmz,eta0);
      }
    }
  }

  if(v>vf_max) {
    alp = 1. - alp;
  }

  return alp;
}
  
