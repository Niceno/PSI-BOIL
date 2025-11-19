#include "vofaxisym.h"

#if 0 /* discarded due to the way elvira works */
/******************************************************************************/
real VOFaxisym::calc_v(const real alpha, const real vma, const real vmb, const real vmc) const {
/***************************************************************************//**
 \brief Solve the point-forward 2D Cartesian problem, i.e. calculate phi(alp,n),
        this is a virtual function wrt 3D VOF calc_v; these two should basically
        produce the same result (2D = special case 3D for Cartesian);
        source: paper by Scardovelli and Zaleski(2000)
*******************************************************************************/

  if(alpha<boil::pico)
    return 0.0;
  if(alpha-1.>-boil::pico)
    return 1.0;

  real v;
  real a;
  real vm1, vm3;
  /* vm2 = vmb = 0! This is 2D xz! */

  a = boil::minr(alpha, 1.0-alpha);
  v = 0.0;

  if(a>0) {
    vm1 = boil::minr(vma,vmc);
    vm3 = 1.0-vm1;

    if(a<vm1) {
      v = a*a/2./vm1/vm3; /* triangle */
    } else {
      v = (a-0.5*vm1)/vm3; /* trapezoid */
    }
  }

  /* problem inversion */
  if(alpha>0.5) {
    v = 1.0-v;
  }

  return v;
}
#endif

/******************************************************************************/
real VOFaxisym::calc_v_axisymmetric(real nnnx, real alp, real eta0, real & K) {
/***************************************************************************//**
 \brief Solve the point-forward axisymmetric problem, i.e. calc phi(alp,n,eta),
        source: derived by me (= Lubomir)
*******************************************************************************/
  if(alp<boil::pico)
    return 0.0;
  if(alp-1.>-boil::pico)
    return 1.0;

  real nnx = nnnx;
  real alpha = alp;

  /* nnx is already normalized by L1 norm */
  real mmx = fabs(nnx);
  real mmz = 1.0-mmx;

  /* complementary problem */
  if(alp>0.5) {
    alpha = 1.0-alpha;
    nnx = -nnx;
  }

  /* Cartesian solution */
  real c = calc_v(alpha,mmx,0.0,mmz);

  /* axisymmetric correction */
  real xi(0.5);
  if(mmx<=mmz) {
    if(nnx>=0.0) {
      if(alpha<mmx) { /* small positive triangle */
        xi = xi_small_pos_triangle(alpha,mmx,mmz);
      } else { /* small positive trapezoid */
        xi = xi_small_pos_trapezoid(alpha,mmx,mmz);
      }
    } else {
      if(alpha<mmx) { /* small negative triangle */
        xi = xi_small_neg_triangle(alpha,mmx,mmz);
      } else { /* small negative trapezoid */
        xi = xi_small_neg_trapezoid(alpha,mmx,mmz);
      }
    }
  } else {
    if(nnx>=0.0) {
      if(alpha<mmz) { /* large positive triangle */
        xi = xi_large_pos_triangle(alpha,mmx,mmz);
      } else { /* large positive trapezoid */
        xi = xi_large_pos_trapezoid(alpha,mmx,mmz);
      }
    } else {
      if(alpha<mmz) { /* large negative triangle */
        xi = xi_large_neg_triangle(alpha,mmx,mmz);
      } else { /* large negative trapezoid */
        xi = xi_large_neg_trapezoid(alpha,mmx,mmz);
      }
    }
  }

  /* correction factor */
  K = (eta0+xi)/(eta0+0.5);

  /* axisymmetric solution */
  real v = c*K;
  //boil::oout<<K<<" "<<v<<" "<<c<<" "<<xi<<boil::endl;

  /* complementary problem */
  if(alp>0.5) {
    v = 1.0-v;
  }

  return v;
}

