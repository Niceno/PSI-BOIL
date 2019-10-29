#include "vof.h"

/******************************************************************************/
void VOF::norm_mixed(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Mixed method in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  real dummy; /* dummy from nalp calculations */

  for_ijk(i,j,k) {

    real nxx, nyy, nzz;
    norm_mixed_kernel(nxx, nyy, nzz, dummy, i,j,k, sca);

    nx[i][j][k] = nxx;
    ny[i][j][k] = nyy;
    nz[i][j][k] = nzz;
  }

  /* boundaries treated using cc */
  norm_cc_near_bnd(sca);

  /* normal vector on boundary plane */
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);

  return;
}

void VOF::norm_mixed_kernel(real & nx_val, real & ny_val, real & nz_val,
                            real & dummy, /* unused */
                            const int i, const int j, const int k,
                            const Scalar & sca) {

  /* normal vector from the CC method; mcomp indicates candidate selected */
  Comp mcomp;
  real nxx_cc, nyy_cc, nzz_cc;
  norm_cc_kernel(nxx_cc, nyy_cc, nzz_cc, dummy, mcomp, i,j,k, sca);

  /* normal vector from the Young method */
  real nxx_young, nyy_young, nzz_young;
  norm_young_kernel(nxx_young, nyy_young, nzz_young, dummy, i,j,k, sca);

  /* selection according to Aulisa (9) */
  bool nyoung = nxx_young*nxx_young+nyy_young*nyy_young+nzz_young*nzz_young < 0.5;
  bool ncc    = nxx_cc*nxx_cc+nyy_cc*nyy_cc+nzz_cc*nzz_cc < 0.5;
  if (nyoung || ncc) {
    nx_val = real(!ncc)*nxx_cc+real(!nyoung)*nxx_young;
    ny_val = real(!ncc)*nyy_cc+real(!nyoung)*nyy_young;
    nz_val = real(!ncc)*nzz_cc+real(!nyoung)*nzz_young;
  } else {
    select_norm_myc(nx_val,ny_val,nz_val, 
                    nxx_cc,nyy_cc,nzz_cc,
                    nxx_young,nyy_young,nzz_young,
                    mcomp);
  }


  return;
}
