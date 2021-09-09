#include "vofaxisym.h"
//#define ONLY_CART
//#define ONLY_CYL

/* axisymmetric curvature is direction-dependent */
real VOFaxisym::calculate_curvature_HF(const arr2D & heights,
                                 const real d1m, const real d1c, const real d1p,
                                 const real d2m, const real d2c, const real d2p,
                                 const bool truedir1, const bool truedir2,
                                 const real mult, const real max_n,
                                 const real xcent) const {

  /* major direction is "i" -> truedir2 = true, heights in the second coord */
  /* major direction is "k" -> truedir1 = true, heights in the first coord */

  real hm, hc, hp, dm, dp;
  bool truedir;
  if(truedir1) {
    //assert(!truedir2); /* sanity check */
    hm = heights[0][1]; 
    hc = heights[1][1]; 
    hp = heights[2][1];
    dm = d1m;
    dp = d1p;
    truedir = truedir1;
  } else {
    //assert(!truedir1); /* sanity check */
    hm = heights[1][0]; 
    hc = heights[1][1]; 
    hp = heights[1][2];
    dm = d2m;
    dp = d2p;
    truedir = truedir2;
  }

  /* Note: Lopez Eq. (6) correction is not used in 2D */
  real h_1c = truedir*(hp-hm)/(dm+dp);
  real h_1u = truedir*(hp-hc)/dp;
  real h_1d = truedir*(hc-hm)/dm;
  real h_11 = truedir*2.*(dp*hm+dm*hp-(dp+dm)*hc)/dp/dm/(dp+dm);

  /* Cartesian contribution */
  real kap_cart = h_11/pow(1.+h_1c*h_1c,1.5);
  real kap_cyl;

  if(truedir1) {
#if 0
    kap_cyl  = 1.0/xcent*h_1c/sqrt(1.+h_1c*h_1c);
#else
    kap_cyl  = h_1u/sqrt(1.+h_1u*h_1u);
    kap_cyl += h_1d/sqrt(1.+h_1d*h_1d);

    kap_cyl /= 2.*xcent;
#endif
  } else {
#if 1
    kap_cyl = -1./xcent * 1./sqrt(1.+h_1c*h_1c);
#else /* this cannot reproduce zero curvature at inflexion */
    kap_cyl  = -1./sqrt(1.+h_1u*h_1u);
    kap_cyl += -1./sqrt(1.+h_1d*h_1d);

    kap_cyl /= 2.*xcent;
#endif
  }

  /* under this convention, bubbles have negative curvature */
  kap_cart *= -mult;
  kap_cyl *= -mult;

#ifdef ONLY_CART
  kap_cyl *= 0.0;
#elif defined ONLY_CYL
  kap_cart *= 0.0;
#endif

  //boil::oout<<hm<<" "<<hc<<" "<<hp<<" | "<<h_1d<<" "<<h_1c<<" "<<h_1u<<" "<<h_11<<" | "<<truedir1<<" | "<<kap_cart<<" "<<kap_cyl<<" "<<kap_cart+kap_cyl<<boil::endl;
  
  return (kap_cyl+kap_cart);
} 
