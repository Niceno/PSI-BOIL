#include "enthalpyfd.h"

/***************************************************************************//**
*  Cell-wise construction of diffusion matrix based on local structure.                              
*******************************************************************************/
void EnthalpyFD::diffmatrix_kernel(const std::array<ConnectType,3> & ctype,
                                   const real cxm, const real cxp,
                                   const std::vector<StencilPoint> & stencil,
                                   real & Am, real & Ac, real & Ap, real & F) {

  /*-----------------------+
  |  fluid - fluid - fluid |
  +-----------------------*/ 
  if(ctype == c_fff) {
    Am = cxm;
    Ap = cxp;
    Ac += cxm + cxp;
  }

  /*---------------------------+
  |  interface - fluid - fluid |
  +---------------------------*/
  else if(ctype == c_iff) {
    //cht.topo->nth_order_first_coefs(coefs,stencil,AccuracyOrder::Second());
    Am = 0.;
    Ap = cxp;
    Ac += cxm + cxp;
    F += cxm*stencil[0].val;
  }

  /*---------------------------+
  |  fluid - fluid - interface |
  +---------------------------*/
  else if(ctype == c_ffi) {
    Am = cxm;
    Ap = 0.;
    Ac += cxm + cxp;
    F += cxp*stencil[2].val;
  }

  /*---------------------------+
  |  fluid - fluid - interface |
  +---------------------------*/
  else if(ctype == c_ifi) {
    Am = 0.;
    Ap = 0.;
    Ac += cxm + cxp;
    F += cxm*stencil[0].val+cxp*stencil[2].val;
  }

  return;
}
