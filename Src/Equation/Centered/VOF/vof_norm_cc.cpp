#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_cc(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Centered Columns method in E.Aulisa,JCP,225(2007),2301-2319
*         Results: nx, ny, nz
*******************************************************************************/

  for_ijk(i,j,k) {
    real nxx, nyy, nzz;
    Comp mcomp;
    norm_cc_kernel(nxx, nyy, nzz, mcomp, i,j,k, sca);

    nx[i][j][k] = nxx;
    ny[i][j][k] = nyy;
    nz[i][j][k] = nzz;
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  //insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
#if 1
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
#endif
  insert_bc_norm_cc(sca);

  /* normalize */
  //for_avijk(sca,i,j,k) {
  //  normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  //}

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);


  return;
}

void VOF::norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val,
                         Comp & mcomp,
                         const int i, const int j, const int k,
                         const Scalar & sca) {

  real nxX, nyX, nzX;
  nxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
  nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
              - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
  nzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
              - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 

  real nxY, nyY, nzY;
  nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
              - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
  nyY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
  nzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));

  real nxZ, nyZ, nzZ;
  nxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
              - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
  nyZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
  nzZ = copysign(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));

  select_norm_cc(nx_val, ny_val, nz_val,
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ,
                 &mcomp);

  return;
}

void VOF::norm_cc_kernel(real & nx_val, real & ny_val, real & nz_val,
                         const int i, const int j, const int k,
                         const Scalar & sca) {

  real nxX, nyX, nzX;
  nxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
  nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
              - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
  nzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
              - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 

  real nxY, nyY, nzY;
  nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
              - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
  nyY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
  nzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));

  real nxZ, nyZ, nzZ;
  nxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
              - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
  nyZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
              - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
  nzZ = copysign(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));

  Comp * mcomp;
  select_norm_cc(nx_val, ny_val, nz_val,
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ,
                 mcomp);

  return;
}
