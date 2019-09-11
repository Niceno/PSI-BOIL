#include "domain.h"

/******************************************************************************/
/* cell surfaces */
real Domain::dSx(const int i, const int j, const int k) const 
  { return dyc(j) * dzc(k); }
real Domain::dSy(const int i, const int j, const int k) const
  { return dxc(i) * dzc(k); }
real Domain::dSz(const int i, const int j, const int k) const
  { return dxc(i) * dyc(j); }

real Domain::dSx_xstag(const int i, const int j, const int k) const
  { return dyc(j) * dzc(k); }
real Domain::dSx_ystag(const int i, const int j, const int k) const
  { return dys(j) * dzc(k); }
real Domain::dSx_zstag(const int i, const int j, const int k) const
  { return dyc(j) * dzb(k); }

real Domain::dSy_xstag(const int i, const int j, const int k) const
  { return dxw(i) * dzc(k); }
real Domain::dSy_ystag(const int i, const int j, const int k) const
  { return dxc(i) * dzc(k); }
real Domain::dSy_zstag(const int i, const int j, const int k) const
  { return dxc(i) * dzb(k); }

real Domain::dSz_xstag(const int i, const int j, const int k) const
  { return dxw(i) * dyc(j); }
real Domain::dSz_ystag(const int i, const int j, const int k) const
  { return dxc(i) * dys(j); }
real Domain::dSz_zstag(const int i, const int j, const int k) const
  { return dxc(i) * dyc(j); }

real Domain::dSx(const Sign sig, const int i, const int j, const int k) const 
  { return dSx(i,j,k); }
real Domain::dSy(const Sign sig, const int i, const int j, const int k) const
  { return dSy(i,j,k); }
real Domain::dSz(const Sign sig, const int i, const int j, const int k) const
  { return dSz(i,j,k); }

real Domain::dSx_xstag(const Sign sig, const int i, const int j, const int k) const
  { return dSx_xstag(i,j,k); }
real Domain::dSx_ystag(const Sign sig, const int i, const int j, const int k) const
  { return dSx_ystag(i,j,k); }
real Domain::dSx_zstag(const Sign sig, const int i, const int j, const int k) const
  { return dSx_zstag(i,j,k); }

real Domain::dSy_xstag(const Sign sig, const int i, const int j, const int k) const
  { return dSy_xstag(i,j,k); }
real Domain::dSy_ystag(const Sign sig, const int i, const int j, const int k) const
  { return dSy_ystag(i,j,k); }
real Domain::dSy_zstag(const Sign sig, const int i, const int j, const int k) const
  { return dSy_zstag(i,j,k); }

real Domain::dSz_xstag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_xstag(i,j,k); }
real Domain::dSz_ystag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_ystag(i,j,k); }
real Domain::dSz_zstag(const Sign sig, const int i, const int j, const int k) const
  { return dSz_zstag(i,j,k); }

/* cell volume */
real Domain::dV(const int i, const int j, const int k) const
  { return dxc(i) * dyc(j) * dzc(k); }

real Domain::dV_xstag(const int i, const int j, const int k) const
  { return dxw(i) * dyc(j) * dzc(k); }
real Domain::dV_ystag(const int i, const int j, const int k) const
  { return dxc(i) * dys(j) * dzc(k); }
real Domain::dV_zstag(const int i, const int j, const int k) const
  { return dxc(i) * dyc(j) * dzb(k); }
