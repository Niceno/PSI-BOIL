#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::advance(const bool anci) {
  advance(phi,anci);
}

void VOFaxisym::advance(Scalar & scp, const bool anci) {
 
  boil::timer.start("vof advance");

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_ijk(i,j,k){
    scp[i][j][k]=scp[i][j][k]+time->dt()*fext[i][j][k];
  }
  scp.bnd_update();
  scp.exchange_all();

  /*-----------------------+
  |  reconstruct geometry  |
  +-----------------------*/
  reconstruct_geometry(scp);

  /* iterate boundary normal vector */

  for_aijk(i,j,k){
    stmp[i][j][k] = scp[i][j][k] * dV(i,j,k);
  }

#if 0
  /* advance in x-direction */
  if(ifull)
    advance_x(scp);
#endif

  /* advance in z-direction */
  if(kfull)
    advance_z(scp);

  /* update phi */
  if (limit_color) {
    for_ijk(i,j,k){
      real phi_tmp = stmp[i][j][k] / dV(i,j,k);
      phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
    }
  } else {
    for_ijk(i,j,k){
      phi[i][j][k] = stmp[i][j][k] / dV(i,j,k);
    }
  }

  phi.bnd_update();
  phi.exchange_all();

  boil::timer.stop("vof advance");

  if(anci)
    ancillary(phi);

  return;
}

