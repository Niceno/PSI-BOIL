#include "vof.h"

/******************************************************************************/
void VOF::advance() {
  advance(phi);
}

void VOF::advance(Scalar & scp) {
 
  boil::timer.start("vof advance");

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_cc(scp);

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_ijk(i,j,k){
    scp[i][j][k]=scp[i][j][k]+time->dt()*fext[i][j][k];
  }
  scp.bnd_update();
  scp.exchange_all();

  for_aijk(i,j,k){
    real phival = scp[i][j][k];
    real fval = fext[i][j][k]*time->dt();

    /* liquid content */
    stmp[i][j][k] = phival * dV(i,j,k);

    /* intermediate phi value */
    /* fext is related to mdot as follows:
     * fext = -m'''/rhol  */
    real denscoef = 1.0-rhol/rhov;
    real denom = std::max(1.0 + denscoef*fval,boil::pico);
    phival /= denom;
    scp[i][j][k] = std::max(0.0,std::min(1.0,phival));
  }

  // advance in x-direction
  advance_x(scp);

  // advance in y-direction
  advance_y(scp);
  
  // advance in z-direction
  advance_z(scp);

  // update phi
  for_ijk(i,j,k){
    real phi_tmp = stmp[i][j][k] / dV(i,j,k);
#if 0
    // limit C
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
#else
    // unlimit C
    phi[i][j][k] = phi_tmp;
#endif

  }
  phi.bnd_update();
  phi.exchange_all();

  boil::timer.stop("vof advance");

  ancillary();

  return;
}

