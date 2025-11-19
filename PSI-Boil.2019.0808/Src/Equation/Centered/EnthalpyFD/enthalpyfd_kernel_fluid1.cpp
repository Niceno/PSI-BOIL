#include "enthalpyfd.h"
#include "def.h"

/***************************************************************************//**
*  Cell-wise construction of diff matrix based on local structure (no solid)
*  - interfacial resistance is considered only in liquid 
*******************************************************************************/
void EnthalpyFD::kernel_fluid1(const std::array<ConnectType,3> & ctype,
                               const real cxm, const real cxp,
                               std::vector<StencilPoint> & stencil,
                               const real resinvm, const real resinvp,
                               real & Am, real & Ac, real & Ap, real & F) {

  /*------------------------+
  |  fluid - fluid - fluid  |
  +------------------------*/ 
  if(ctype == c_fff) {
    Am = cxm;
    Ap = cxp;
    Ac += cxm + cxp;
  }

  /*----------------------------+
  |  interface - fluid - fluid  |
  +----------------------------*/
  else if(ctype == c_iff) {
    Am = 0.;
    Ap = cxp;
    Ac += cxm + cxp;

    /* do we use interfacial heat transfer resistance? */
    if(boil::realistic(resinvm)) {
      /* correct for upwind scheme */
      stencil[2].pos -= stencil[0].pos;
      stencil[1].pos -= stencil[0].pos;
      stencil[0].pos -= stencil[0].pos;

      /* evaluate coefs */
      std::vector<real> coefs;
#ifdef USE_FIRST_ORDER_INTRESIST
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::First());

      /* correct matrix */
      //Ap += 0.;
#else
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());

      /* correct matrix */
      Ap +=  cxm*coefs[2]/(resinvm-coefs[0]);
#endif
      Ac += -cxm*coefs[1]/(resinvm-coefs[0]);
      F += cxm*stencil[0].val*resinvm/(resinvm-coefs[0]);  
    } else {
      F += cxm*stencil[0].val;
    }
  }

  /*----------------------------+
  |  fluid - fluid - interface  |
  +----------------------------*/
  else if(ctype == c_ffi) {
    Am = cxm;
    Ap = 0.;
    Ac += cxm + cxp;

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
#ifdef USE_FIRST_ORDER_INTRESIST
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::First());

      /* correct matrix */
      //Am +=  0.;
#else
      cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());

      /* correct matrix */
      Am +=  cxp*coefs[2]/(resinvp-coefs[0]);
#endif
      Ac += -cxp*coefs[1]/(resinvp-coefs[0]);
      F += cxp*stencil[0].val*resinvp/(resinvp-coefs[0]);
    } else {
      /* here the stencil is not reversed */
      F += cxp*stencil[2].val;
    }
  }

  /* this one is always first-order discretisation */

  /*--------------------------------+
  |  interface - fluid - interface  |
  +--------------------------------*/
  else if(ctype == c_ifi) {
    Am = 0.;
    Ap = 0.;
    Ac += cxm + cxp;

    if(boil::realistic(resinvp)) {
      /* sanity check */
      assert(boil::realistic(resinvm));

      /* we default to first order for simplicity */
      std::vector<StencilPoint> sm, sp;
      sm.push_back(stencil[0]);
      sm.push_back(stencil[1]);

      sp.push_back(stencil[2]);
      sp.push_back(stencil[1]);

      /* reverse the stencil */
      sp[0].pos *= -1.;

      /* correct for upwind scheme */
      sm[1].pos -= sm[0].pos;
      sm[0].pos -= sm[0].pos;
      sp[1].pos -= sp[0].pos;
      sp[0].pos -= sp[0].pos;

      /* evaluate coefs */
      std::vector<real> cms, cps;
      cht.topo->nth_order_first_coefs(cms,sm,AccuracyOrder::First());
      cht.topo->nth_order_first_coefs(cps,sp,AccuracyOrder::First());

      /* correct matrix */
      Ac += -cxm*cms[1]/(resinvm-cms[0])
            -cxp*cps[1]/(resinvp-cps[0]);
      F += cxm*sm[0].val*resinvm/(resinvm-cms[0])
         + cxp*sp[0].val*resinvp/(resinvp-cps[0]);
    } else {
      F += cxm*stencil[0].val+cxp*stencil[2].val;
    }
  }

  else {
    boil::aout<<"EFD::diffmatrix: Incorrect kernel requested! Exiting."
              <<boil::endl;
    exit(0);
  }

  return;
}
