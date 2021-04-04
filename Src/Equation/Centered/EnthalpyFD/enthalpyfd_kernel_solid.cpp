#include "enthalpyfd.h"

/***************************************************************************//**
*  Cell-wise construction of diff matrix based on local structure (in solid)
*******************************************************************************/
void EnthalpyFD::kernel_solid(const std::array<ConnectType,3> & ctype,
                              const real cxm, const real cxp,
                              const real pm, const real pp,
                              const std::array<real,3> resistvals,
                              const real dwsrcm, const real dwsrcp,
                              real & Am, real & Ac, real & Ap, real & F) {

  /*------------------------+
  |  solid - solid - solid  |
  +------------------------*/ 
  if(ctype == c_sss) {
    Am = cxm;
    Ap = cxp;
    Ac += cxm + cxp;
  }

  /*----------------------------------+
  |  interface/fluid - solid - solid  |
  +----------------------------------*/
  else if(ctype == c_iss || ctype == c_fss) {
    real fact = resistance_multiplier(resistvals[1],resistvals[0]);
    Ap = cxp;
    Ac += fact*cxm + cxp;
    if(ctype == c_fss) {
      Am = fact*cxm;
    } else {
      Am = 0.;
      F += fact*cxm*pm;
    }
    F += dwsrcm*fact*cxm*resistvals[0];
  }

  /*----------------------------------+
  |  solid - solid - interface/fluid  |
  +----------------------------------*/
  else if(ctype == c_ssi || ctype == c_ssf) {
    real fact = resistance_multiplier(resistvals[1],resistvals[2]);
    Am = cxm;
    Ac += cxm + fact*cxp;
    if(ctype == c_ssf) {
      Ap = fact*cxp;
    } else {
      Ap = 0.;
      F += fact*cxp*pp;
    }
    F += dwsrcp*fact*cxp*resistvals[2];
  }

  else {
    boil::aout<<"EFD::diffmatrix: Incorrect kernel requested! Exiting."
              <<boil::endl;
    exit(0);
  }

  return;
}
