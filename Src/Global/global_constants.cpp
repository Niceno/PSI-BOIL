#include "global_constants.h"

namespace boil {
  /* big */
  real yotta(1.0e+24);
  real zetta(1.0e+21);
  real exa  (1.0e+18);
  real peta (1.0e+15);
  real tera (1.0e+12); 
  real giga (1.0e+09); 
  real mega (1.0e+06);
  real kilo (1.0e+03);
  /* small */
  real milli(1.0e-03);
  real micro(1.0e-06);
  real nano (1.0e-09); 
  real pico (1.0e-12); 
  real femto(1.0e-15);
  real atto (1.0e-18);
  real zepto(1.0e-21);
  real yocto(1.0e-24);

  real unreal = yotta;

  real pi(acos(-1.0));

  /* gravitational constant for earth.  
     should be changed if martians start to use psi-boil */
  real g(9.81);
}
