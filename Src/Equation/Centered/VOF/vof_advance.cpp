#include "vof.h"

/******************************************************************************/
void VOF::advance(const bool anci) {
  advance(phi,anci);
}

void VOF::advance(Scalar & scp, const bool anci) {
 
  boil::timer.start("vof advance");

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_ijk(i,j,k){
    scp[i][j][k]=scp[i][j][k]+time->dt()*fext[i][j][k];
  }
  scp.bnd_update();
  scp.exchange_all();

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_mixed(scp);
  //norm_cc(scp);
  //norm_elvira(scp);

  for_aijk(i,j,k){
    stmp[i][j][k] = scp[i][j][k] * dV(i,j,k);
  }

  // advance in x-direction
  if(ifull)
    advance_x(scp);

  // advance in y-direction
  if(jfull)
    advance_y(scp);
  
  // advance in z-direction
  if(kfull)
    advance_z(scp);

  // update phi
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
    ancillary();

  return;
}

