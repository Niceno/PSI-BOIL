#include "enthalpyfd.h"

/***************************************************************************//**
*  Cell-wise construction of diff matrix based on local structure (fluid at sol)
*  - interfacial resistance is considered only in liquid 
*******************************************************************************/
void EnthalpyFD::kernel_fluid2(const std::array<ConnectType,3> & ctype,
                               const real cxm, const real cxp,
                               std::vector<StencilPoint> & stencil,
                               const real resinvm, const real resinvp,
                               const std::array<real,3> resistvals,
                               const real dwsrcm, const real dwsrcp,
                               real & Am, real & Ac, real & Ap, real & F) {

  /*------------------------+
  |  solid - fluid - fluid  |
  +------------------------*/
  if(ctype == c_sff) {
    real fact = resistance_multiplier(resistvals[1],resistvals[0]);
    Am = fact*cxm;
    Ap = cxp;
    Ac += fact*cxm + cxp;
    F += dwsrcm*fact*cxm*resistvals[0];
  }

  /*------------------------+
  |  fluid - fluid - solid  |
  +------------------------*/
  else if(ctype == c_ffs) {
    real fact = resistance_multiplier(resistvals[1],resistvals[2]);
    Am = cxm;
    Ap = fact*cxp;
    Ac += cxm + fact*cxp;
    F += dwsrcp*fact*cxp*resistvals[2];
  }

  /*----------------------------+
  |  interface - fluid - solid  |
  +----------------------------*/
  else if(ctype == c_ifs) {
    real fact = resistance_multiplier(resistvals[1],resistvals[2]);
    Am = 0.;
    Ap = fact*cxp;
    Ac += cxm + fact*cxp;
    F += dwsrcp*fact*cxp*resistvals[2];

    /* do we use interfacial heat transfer resistance? */
    if(boil::realistic(resinvm)) {
      /* correct for upwind scheme */
      stencil[2].pos -= stencil[0].pos;
      stencil[1].pos -= stencil[0].pos;
      stencil[0].pos -= stencil[0].pos;

      /* evaluate coefs */
      std::vector<real> coefs;
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());

      /* correct matrix */
      real deno = resinvm-coefs[0];
      real mult = coefs[2]/deno;
      real fact2 = resistance_multiplier(resistvals[2],resistvals[1]);
      Ap +=  cxm * mult * fact;
      Ac += -cxm * (coefs[1]/deno + mult*fact2);
      F += cxm*stencil[0].val*resinvm/deno
         + dwsrcp * cxm * mult * fact * resistvals[2];

    } else {
      F += cxm*stencil[0].val;
    }
  }

  /*----------------------------+
  |  solid - fluid - interface  |
  +----------------------------*/
  else if(ctype == c_sfi) {
    real fact = resistance_multiplier(resistvals[1],resistvals[0]);
    Am = fact*cxm;
    Ap = 0.;
    Ac += fact*cxm + cxp;
    F += dwsrcm*fact*cxm*resistvals[0];

    /* do we use interfacial heat transfer resistance? */
    if(boil::realistic(resinvp)) {
      /* reverse the stencil */
      std::reverse(stencil.begin(),stencil.end());
      stencil[2].pos *= -1.;
      stencil[0].pos *= -1.;

      /* correct for upwind scheme */
      stencil[2].pos -= stencil[0].pos;
      stencil[1].pos -= stencil[0].pos;
      stencil[0].pos -= stencil[0].pos;

      /* evaluate coefs */
      std::vector<real> coefs;
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());

      /* correct matrix */
      real deno = resinvp-coefs[0];
      real mult = coefs[2]/deno;
      real fact2 = resistance_multiplier(resistvals[0],resistvals[1]);
      Am +=  cxp * mult * fact;
      Ac += -cxp * (coefs[1]/deno + mult*fact2);
      F += cxp*stencil[0].val*resinvp/deno
         + dwsrcm * cxp * mult * fact * resistvals[2]; 

    } else {
      /* here the stencil is not reversed */
      F += cxp*stencil[2].val;
    }
  }

  else {
    boil::aout<<"EFD::diffmatrix: Incorrect kernel requested! Exiting."
              <<boil::endl;
    exit(0);
  }

  return;
}
