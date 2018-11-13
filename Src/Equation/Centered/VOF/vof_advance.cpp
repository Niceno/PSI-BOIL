#include "vof.h"

/******************************************************************************/
void VOF::advance() {

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_aijk(i,j,k){   //must be aijk for insert boundary
    phi[i][j][k]=phi[i][j][k]+time->dt()*fext[i][j][k];
  }
  phi.bnd_update();
  phi.exchange_all();

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  //gradphic(phi);
  norm_cc(phi);

  for_aijk(i,j,k){
    stmp[i][j][k] = phi[i][j][k] * dV(i,j,k);
  }

  // advance in x-direction
  advance_x();

  // advance in y-direction
  advance_y();
  
  // advance in z-direction
  advance_z();

  // update phi
  for_ijk(i,j,k){
    real phi_tmp = stmp[i][j][k] / dV(i,j,k);
#if 0
    // limit C
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
    if(phi_tmp>1.0+boil::pico || phi_tmp< -boil::pico){
      std::cout.setf(std::ios_base::scientific);
      std::cout<<"limit phi "<<phi_tmp<<" i "<<i<<" j "<<j<<" k "<<k<<"\n";
      std::cout.unsetf(std::ios_base::floatfield);
    }
#else
    // unlimit C
    phi[i][j][k] = phi_tmp;
#endif

  }
  phi.bnd_update();
  phi.exchange_all();

  return;
}

