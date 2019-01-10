#include "vof.h"

/******************************************************************************/
void VOF::advance() {

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_aijk(i,j,k){   //must be aijk for insert boundary
   #if 1
    /* fext cut */
    real ftot = fext[i][j][k];
    real fval = fext_cut(i,j,k,ftot);
    phi[i][j][k]=phi[i][j][k]+time->dt()*fval;
    stmp2[i][j][k] = ftot-fval;
   #else
    phi[i][j][k]=phi[i][j][k]+time->dt()*fext[i][j][k];
   #endif

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

#if 1
  adens.exchange();
  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();
#endif
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

  #if 1
    /* fext cut II */
    real ftot = stmp2[i][j][k];
    real fval = fext_cut(i,j,k,ftot);
    phi[i][j][k]=phi[i][j][k]+time->dt()*fval;
  #endif
#endif
  }
  phi.bnd_update();
  phi.exchange_all();

#if 1
  /* calculate alpha in cells */
  extract_alpha();

  /*  calculate free surface position */
  cal_fs3();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate area */
  cal_adens();

  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();
#endif

#if 0
  curv_HF();
#endif
  return;
}

