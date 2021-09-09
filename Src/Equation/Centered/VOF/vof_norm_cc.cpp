#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_cc(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Centered Columns method in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  real dummy; /* dummy from alpha calculations */

  for_ijk(i,j,k) {
    real nxx, nyy, nzz;
    Comp mcomp;
    norm_cc_kernel(nxx, nyy, nzz, dummy, mcomp, i,j,k, sca);

    nx[i][j][k] = nxx;
    ny[i][j][k] = nyy;
    nz[i][j][k] = nzz;
  }

  /* normal vector on boundary plane */
  norm_cc_near_bnd(sca);
#if 1
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
#endif

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);


  return;
}

void VOF::norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val,
                         real & dummy, /* unused */
                         Comp & mcomp,
                         const int i, const int j, const int k,
                         const Scalar & sca) {

  real nxX, nyX, nzX;
  nxX = signum(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
  nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
              - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
  nzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
              - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 

  real nxY, nyY, nzY;
  nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
              - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
  nyY = signum(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
  nzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));

  real nxZ, nyZ, nzZ;
  nxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
              - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
  nyZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
  nzZ = signum(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));

  select_norm_cc(nx_val, ny_val, nz_val,
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ,
                 &mcomp);

  return;
}

void VOF::norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val,
                         real & dummy, /* unused */
                         const int i, const int j, const int k,
                         const Scalar & sca) {

  real nxX, nyX, nzX;
  nxX = signum(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
  nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
              - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
  nzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
              - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 

  real nxY, nyY, nzY;
  nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
              - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
  nyY = signum(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
  nzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));

  real nxZ, nyZ, nzZ;
  nxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
              - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
  nyZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
  nzZ = signum(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));

  Comp mcomp;
  select_norm_cc(nx_val, ny_val, nz_val,
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ,
                 &mcomp);

  return;
}
