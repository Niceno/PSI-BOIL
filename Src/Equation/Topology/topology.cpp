#include "topology.h"

Topology::Topology(Scalar * VF, Scalar * CLR,
                   Scalar * NX, Scalar * NY, Scalar * NZ, 
                   Scalar * ADENS, Vector * FS, ScalarInt * IFLAG,
                   const real CLRSURF) :
  vf(VF),
  clr(CLR),
  nx(NX),
  ny(NY),
  nz(NZ),
  fs(FS),
  adens(ADENS),
  iflag(IFLAG),
  stmp(*CLR->domain()),
  delta(*CLR->domain()),
  stmp2(*CLR->domain()),
  iflagold(*IFLAG->domain()),
  clrold(*CLR->domain()), 
  vfold(*VF->domain()),
  fsold(*FS->domain()),
  bflag_struct(*CLR),
  clrsurf(CLRSURF) {

  stmp  = NX->shape();
  stmp2 = NX->shape();
  delta = NX->shape();

  iflagold = IFLAG->shape();
  clrold = CLR->shape();
  vfold = VF->shape();
  for_m(m) {
    fsold(m) = (*FS)(m).shape();
  }

  dxmin = domain()->dxyz_min();
  mmax_ext = 100;
  tol_ext = 1e-7; 
  close_to_cc = 1.0e-2;
}
