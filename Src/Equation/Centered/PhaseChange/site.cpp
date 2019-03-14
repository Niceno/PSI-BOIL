#include "site.h"

/******************************************************************************/
Site::Site (const real x, const real y, const real z,
            const real t, const real zpl) {
  xsite=x;
  ysite=y;
  if (z<boil::nano) {
    //boil::oout<<"site: z is moved to "<<boil::pico<<"\n";
    zsite=boil::nano;
  } else {
    zsite=z;
  }
  set_active_tpr(t);
  set_zplant(zpl);
  tseed=-boil::exa;
  tcutneck=-boil::exa;
  bseed_prev = false;
  bneck_prev = false;
  req_seed   =true;  // request from site for replant
  allow_seed =true;  // answer to site for replant
  set_active(false); // default value. It is modified in phasechange_micro
  bqsink = false;

#if 0
  real h = z + rs;                             // height of bubble seed
  real r_base = sqrt(h*(2.0*rs-h));            // raduis of bubble-base
  vs = boil::pi/6.0*h*(3.0*r_base*r_base+h*h); // volume of seed
  ab = boil::pi*r_base*r_base;                 // area of bubble-base
#endif
}

Site::~Site () {
}
