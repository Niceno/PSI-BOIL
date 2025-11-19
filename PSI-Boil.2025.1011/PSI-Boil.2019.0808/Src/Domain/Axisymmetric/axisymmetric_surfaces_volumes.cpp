#include "axisymmetric.h"

/***************************************************************************//**
* A derived class from the Domain class
* - x-direction is the radial direction
* - z-direction is the axial direction
* - y-direction is a dummy direction -> dSy is 0!
* - a wedge with angle 1 rad is assumed in the calculations
* - fabs to prevent negative values due to symmetry at origin
*******************************************************************************/
/* cell surfaces */
real Axisymmetric::dSx(const int i, const int j, const int k) const 
  { return fabs(xc(i)) * dzc(k); }  /* area of cylindrical wedge in radial dir */
real Axisymmetric::dSy(const int i, const int j, const int k) const
  { return 0.0; } /* dummy direction */
real Axisymmetric::dSz(const int i, const int j, const int k) const
  { return 0.5*fabs(xn(i+1)*xn(i+1)-xn(i)*xn(i)); } // area of 1 rad
  /* truncated circular sector, angle of 1 rad */

real Axisymmetric::dSx_xstag(const int i, const int j, const int k) const
  //{ return (xc(i)-0.5*dxw(i)) * dzc(k); } /* doesn't work @ origin! */
  { return fabs(xn(i)) * dzc(k); } 
real Axisymmetric::dSx_ystag(const int i, const int j, const int k) const
  { return dSx(i,j,k); } /* pseudo-staggering */
real Axisymmetric::dSx_zstag(const int i, const int j, const int k) const
  { return fabs(xc(i)) * dzb(k); }

real Axisymmetric::dSy_xstag(const int i, const int j, const int k) const
  { return 0.0; }
real Axisymmetric::dSy_ystag(const int i, const int j, const int k) const
  { return 0.0; }
real Axisymmetric::dSy_zstag(const int i, const int j, const int k) const
  { return 0.0; }

real Axisymmetric::dSz_xstag(const int i, const int j, const int k) const
  /* this is so complicated to safeguard @ origin against 0 area */
  // area of 1 rad
  { return 0.5*(fabs(xc(i)*xc(i)-xn(i)*xn(i))+fabs(xn(i)*xn(i)-xc(i-1)*xc(i-1))); }
real Axisymmetric::dSz_ystag(const int i, const int j, const int k) const
  { return dSz(i,j,k); } /* pseudo-staggering */
real Axisymmetric::dSz_zstag(const int i, const int j, const int k) const
  { return dSz(i,j,k); } /* stag in z doesn't change area */

real Axisymmetric::dSx(const Sign sig, const int i, const int j, const int k) const 
  { return sig == Sign::neg() ? fabs(xn(i))* dzc(k) : fabs(xn(i+1)) * dzc(k); }
real Axisymmetric::dSy(const Sign sig, const int i, const int j, const int k) const
  { return 0.0; }
real Axisymmetric::dSz(const Sign sig, const int i, const int j, const int k) const
  { return dSz(i,j,k); } /* z area doesn't change */

real Axisymmetric::dSx_xstag(const Sign sig, const int i, const int j, const int k) const
  { return sig == Sign::neg() ? fabs(xc(i-1)) * dzc(k) : fabs(xc(i)) * dzc(k); }
real Axisymmetric::dSx_ystag(const Sign sig, const int i, const int j, const int k) const
  { return dSx(sig,i,j,k); } /* pseudo-staggering */
real Axisymmetric::dSx_zstag(const Sign sig, const int i, const int j, const int k) const
  { return sig == Sign::neg() ? fabs(xn(i)) * dzb(k) : fabs(xn(i+1)) * dzb(k); }

real Axisymmetric::dSy_xstag(const Sign sig, const int i, const int j, const int k) const
  { return 0.0; }
real Axisymmetric::dSy_ystag(const Sign sig, const int i, const int j, const int k) const
  { return 0.0; }
real Axisymmetric::dSy_zstag(const Sign sig, const int i, const int j, const int k) const
  { return 0.0; }

  /* sign does not affect z area */
real Axisymmetric::dSz_xstag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_xstag(i,j,k); }
real Axisymmetric::dSz_ystag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_ystag(i,j,k); }
real Axisymmetric::dSz_zstag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_zstag(i,j,k); }

/* cell volume */
real Axisymmetric::dV(const int i, const int j, const int k) const
  { return dSz(i,j,k) * dzc(k); }

real Axisymmetric::dV_xstag(const int i, const int j, const int k) const
  { return dSz_xstag(i,j,k) * dzc(k); }
real Axisymmetric::dV_ystag(const int i, const int j, const int k) const
  { return dSz_ystag(i,j,k) * dzc(k); }
real Axisymmetric::dV_zstag(const int i, const int j, const int k) const
  { return dSz_zstag(i,j,k) * dzb(k); }

/* area value for Cartesian projection */
real Axisymmetric::dSy_cartesian(const int i, const int j, const int k) const
  { return dxc(i)*dzc(k); }

