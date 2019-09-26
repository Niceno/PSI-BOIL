#include "vof.h"

/******************************************************************************/
void VOF::norm_mixed(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Mixed method in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  for_ijk(i,j,k) {

    real nxx, nyy, nzz;
    norm_mixed_kernel(nxx, nyy, nzz, i,j,k, sca);

    nx[i][j][k] = nxx;
    ny[i][j][k] = nyy;
    nz[i][j][k] = nzz;
  }

  /* normal vector on boundary plane */
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
  /* boundaries treated using cc */
  insert_bc_norm_cc(sca);

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);

  return;
}

void VOF::norm_mixed_kernel(real & nx_val, real & ny_val, real & nz_val,
                            const int i, const int j, const int k,
                            const Scalar & sca) {

  /* normal vector from the CC method; mcomp indicates candidate selected */
  Comp mcomp;
  real nxx_cc, nyy_cc, nzz_cc;
  norm_cc_kernel(nxx_cc, nyy_cc, nzz_cc, mcomp, i,j,k, sca);

  /* normal vector from the Young method */
  real nxx_young, nyy_young, nzz_young;
  norm_young_kernel(nxx_young, nyy_young, nzz_young, i,j,k, sca);

  /* selection according to Aulisa (9) */
  real nxx, nyy, nzz;
  select_norm_myc(nxx,nyy,nzz, 
                  nxx_cc,nyy_cc,nzz_cc,
                  nxx_young,nyy_young,nzz_young,
                  mcomp);

  nx_val = nxx;
  ny_val = nyy;
  nz_val = nzz;

  return;
}
