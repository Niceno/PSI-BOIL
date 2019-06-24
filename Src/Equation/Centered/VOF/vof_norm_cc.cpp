#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_cc(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Centered Columns method in E.Aulisa,JCP,225(2007),2301-2319
*         Resluts: nx, ny, nz
*******************************************************************************/

  for_ijk(i,j,k) {
#if 1
    real nxX, nyX, nzX;
    nxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
    nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
                - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
    nzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
                - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 
    normalize(nxX,nyX,nzX);

    real nxY, nyY, nzY;
    nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
                - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
    nyY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
    nzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));
    normalize(nxY,nyY,nzY);

    real nxZ, nyZ, nzZ;
    nxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
                - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
    nyZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
    nzZ = copysign(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));
    normalize(nxZ,nyZ,nzZ);

    if (fabs(nxX)<fabs(nyY)) {
      if (fabs(nyY)<fabs(nzZ)) {
        nx[i][j][k]=nxZ;
        ny[i][j][k]=nyZ;
        nz[i][j][k]=nzZ;
      } else {
        nx[i][j][k]=nxY;
        ny[i][j][k]=nyY;
        nz[i][j][k]=nzY;
      }
    } else {
      if (fabs(nxX)<fabs(nzZ)) {
        nx[i][j][k]=nxZ;
        ny[i][j][k]=nyZ;
        nz[i][j][k]=nzZ;
      } else {
        nx[i][j][k]=nxX;
        ny[i][j][k]=nyX;
        nz[i][j][k]=nzX;
      }
    }
#else
        nx[i][j][k]=phi.xc(i);
        ny[i][j][k]=phi.yc(j);
        nz[i][j][k]=0.0;

        normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
#endif
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  //insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
  insert_bc_norm_cc(sca);
  
  /* normalize */
  //for_avijk(sca,i,j,k) {
  //  normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  //}
 
#if 0
  for_avijk(sca,i,j,k)
    if(k==boil::BW && ny[i][j][k] != 0.0)
      boil::oout<<i<<" "<<j<<" "<<sca[i][j][k]<<" "<<nx[i][j][k]<<" "<<ny[i][j][k]<<" "<<nz[i][j][k]<<boil::endl;
  exit(0);
#endif

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);


  return;
}
