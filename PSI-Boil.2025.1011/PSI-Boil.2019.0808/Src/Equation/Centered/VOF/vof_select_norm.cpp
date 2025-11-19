#include "vof.h"

/******************************************************************************/
void VOF::select_norm_cc(real & nx_val, real & ny_val, real & nz_val,
                         real & nxX, real & nyX, real & nzX,
                         real & nxY, real & nyY, real & nzY,
                         real & nxZ, real & nyZ, real & nzZ,
                         Comp * mcomp) {
/***************************************************************************//**
*  \brief Select normal vector according to CC criterion
*  in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  real nxX_l1, nyX_l1, nzX_l1;
  normalize_l1(nxX_l1, nyX_l1, nzX_l1 , nxX,nyX,nzX);
  normalize(nxX,nyX,nzX);

  real nxY_l1, nyY_l1, nzY_l1;
  normalize_l1(nxY_l1, nyY_l1, nzY_l1 , nxY,nyY,nzY);
  normalize(nxY,nyY,nzY);

  real nxZ_l1, nyZ_l1, nzZ_l1;
  normalize_l1(nxZ_l1, nyZ_l1, nzZ_l1 , nxZ,nyZ,nzZ);
  normalize(nxZ,nyZ,nzZ);

#if 0 /* incorrect selection */
  if (fabs(nxX)<fabs(nyY)) {
    if (fabs(nyY)<fabs(nzZ)) {
      nx_val=nxZ;
      ny_val=nyZ;
      nz_val=nzZ;
    } else {
      nx_val=nxY;
      ny_val=nyY;
      nz_val=nzY;
    }
  } else {
    if (fabs(nxX)<fabs(nzZ)) {
      nx_val=nxZ;
      ny_val=nyZ;
      nz_val=nzZ;
    } else {
      nx_val=nxX;
      ny_val=nyX;
      nz_val=nzX;
    }
  }
#else
  if (fabs(nxX_l1)<fabs(nyY_l1)) {
    if (fabs(nyY_l1)<fabs(nzZ_l1)) {
      nx_val=nxZ;
      ny_val=nyZ;
      nz_val=nzZ;
      *mcomp = Comp::k();
    } else {
      nx_val=nxY;
      ny_val=nyY;
      nz_val=nzY;
      *mcomp = Comp::j();
    }
  } else {
    if (fabs(nxX_l1)<fabs(nzZ_l1)) {
      nx_val=nxZ;
      ny_val=nyZ;
      nz_val=nzZ;
      *mcomp = Comp::k();
    } else {
      nx_val=nxX;
      ny_val=nyX;
      nz_val=nzX;
      *mcomp = Comp::i();
    }
  }
#endif

  return;
}

/******************************************************************************/
void VOF::select_norm_myc(real & nx_val, real & ny_val, real & nz_val,
                          const real & nx_cc, const real & ny_cc, const real & nz_cc,
                          const real & nx_young, const real & ny_young, const real & nz_young,
                          const Comp & mcomp) {
/***************************************************************************//**
*  \brief Select normal vector according to MYC criterion
*  in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  /* selection is done based on the L1 norm */
  real nx_cc_l1, ny_cc_l1, nz_cc_l1;
  normalize_l1(nx_cc_l1, ny_cc_l1, nz_cc_l1 , nx_cc,ny_cc,nz_cc);

  real nx_young_l1, ny_young_l1, nz_young_l1;
  normalize_l1(nx_young_l1, ny_young_l1, nz_young_l1 , nx_young,ny_young,nz_young);

  if       (mcomp==Comp::i()) {

    if(fabs(nx_cc_l1)<fabs(nx_young_l1)) {
      nx_val=nx_cc;
      ny_val=ny_cc;
      nz_val=nz_cc;
    } else {
      nx_val=nx_young;
      ny_val=ny_young;
      nz_val=nz_young;
    }

  } else if(mcomp==Comp::j()) {

    if(fabs(ny_cc_l1)<fabs(ny_young_l1)) {
      nx_val=nx_cc;
      ny_val=ny_cc;
      nz_val=nz_cc;
    } else {
      nx_val=nx_young;
      ny_val=ny_young;
      nz_val=nz_young;
    }

  } else {

    if(fabs(nz_cc_l1)<fabs(nz_young_l1)) {
      nx_val=nx_cc;
      ny_val=ny_cc;
      nz_val=nz_cc;
    } else {
      nx_val=nx_young;
      ny_val=ny_young;
      nz_val=nz_young;
    }

  }

  return;
}
