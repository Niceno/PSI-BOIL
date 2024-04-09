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
#if 1
  #if 1
  norm_mixed(scp);
  #else
  norm_cc(scp);
  #endif
#else
  norm_elvira(scp);
#endif
  bdnorm(scp);

#if 0
  boil::plot->plot(scp,nx,ny,nz, "vof-nx-ny-nz", time->current_step());
  exit(0);
#endif

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
    int ierr=0;
    for_ijk(i,j,k){
      phi[i][j][k] = stmp[i][j][k] / dV(i,j,k);
      if (phi[i][j][k]<-1.0||phi[i][j][k]>2.0) {
        boil::aout<<"Error!!! Too small or too large phi in vof_advance. phi= "
                  <<phi[i][j][k]<<"\n";
        boil::aout<<"proc= "<<boil::cart.iam()
                  <<" (i,j,k)= "<<i<<" "<<j<<" "<<k<<"\n";
        boil::aout<<"(x,y,z)= "<<phi.xc(i)<<" "<<phi.yc(j)<<" "<<phi.zc(k)<<"\n";
        ierr=1;
      }
    }
    boil::cart.sum_int(&ierr);
    if (ierr>=1) {
      exit(0);
    }
  }

  phi.bnd_update();
  phi.exchange_all();

  boil::timer.stop("vof advance");

  if(anci)
    ancillary(phi);

  return;
}

