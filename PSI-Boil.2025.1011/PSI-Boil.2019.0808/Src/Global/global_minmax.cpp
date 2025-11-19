#include "global_precision.h"
#include "global_minmax.h"

namespace boil {

/******************************************************************************/
real maxr(const real & a, const real & b) {
  if(a > b) 
    return a;
  else 
    return b;
}

/******************************************************************************/
real maxr(const real & a, const real & b, const real & c) {
  return maxr(a, maxr(b,c));
}

/******************************************************************************/
real minr(const real & a, const real & b) {
  if(a < b) 
    return a;
  else 
    return b;
}

/******************************************************************************/
real minr(const real & a, const real & b, const real & c) {
  return minr(a, minr(b,c));
}

/******************************************************************************/
int maxi(const int & a, const int & b) {
  if(a > b) 
    return a;
  else 
    return b;
}

/******************************************************************************/
int maxi(const int & a, const int & b, const int & c) {
  return maxi(a, maxi(b,c));
}

/******************************************************************************/
int mini(const int & a, const int & b) {
  if(a < b) 
    return a;
  else 
    return b;
}

/******************************************************************************/
int mini(const int & a, const int & b, const int & c) {
  return mini(a, mini(b,c));
}

} /* boil */
