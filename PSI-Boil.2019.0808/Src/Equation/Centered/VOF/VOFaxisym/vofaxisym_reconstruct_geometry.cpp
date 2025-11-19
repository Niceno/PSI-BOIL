#include "vofaxisym.h"

void VOFaxisym::reconstruct_geometry() {
  reconstruct_geometry(phi);
  return;
}

/******************************************************************************/
void VOFaxisym::reconstruct_geometry(Scalar & scp) {
/***************************************************************************//**
 \brief Reconstruct the geometry of the interface, i.e. calculate the normal
        vector, line constant alpha and color function. The reconstruction pro-
        cess is iterative.
    line: vm1*x + vm2*z = alpha
    output: nx, ny = 0.0, nz, nalpha, c
*******************************************************************************/

  real linferr(boil::unreal);
  int iter(0);

  /* direct iteration */

  /* set c0 = phi */
  for_avijk(clr,i,j,k) {
    clr[i][j][k] = std::max(0.0,std::min(1.0,scp[i][j][k]));
  }

  /* true = extract alpha */
  norm(clr,norm_method_advance,true);

  /* iterate boundary normal vector */
  bdnorm(clr);

  /* doing the reconstruction in the whole domain is unstable; therefore,
     only a band around the presumed interface is considered */
  const int nlayer = 2; //4
  set_reconstruction_flag(scp,tempflag,nlayer);

#if 1
  /* iteration procedure */
  do {
    backward_axisymmetric(scp,nalpha);
    //forward_cartesian(axistmp);
    forward(axistmp);

    /* correction */
    for_avijk(axistmp,i,j,k) {
      if(abs(tempflag[i][j][k])>nlayer) {
        axistmp[i][j][k] = std::max(0.0,std::min(1.0,scp[i][j][k]));
      }
    }

    /* color to vf contains calculations of normal vector and update at walls */
    color_to_vf(axistmp,stmp2,true,true,true);
    /* calculate Linf error norm of reconstruction in phi-space */
  #if 0
    int imax,jmax,kmax;
    bool elvibool;
  #endif
    linferr = linf_scalar_error(stmp2,scp);
                                //,imax,jmax,kmax,elvibool);

    iter++;
  #if 0
    boil::oout<<"VOFaxisym::reconstruct_geometry: "<<time->current_time()
              <<" "<<iter<<" "<<linferr
              <<" | "<<" "<<imax<<" "<<kmax<<" "
              <<" | "<<axistmp[imax][jmax][kmax]<<" "<<stmp2[imax][jmax][kmax]<<" "<<scp[imax][jmax][kmax]<<" "<<nalpha[imax][jmax][kmax]/(fabs(nx[imax][jmax][kmax])+fabs(nz[imax][jmax][kmax])+boil::atto)<<" "<<nx[imax][jmax][kmax]<<" "<<nz[imax][jmax][kmax]<<" "<<elvibool
              <<boil::endl;
    boil::plot->plot(scp,axistmp,nx,nz,nalpha,"phi-clr-nx-nz-nalp", time->current_step());
    if(time->current_step()==1)exit(0);
  #endif

    for_avijk(clr,i,j,k) {
      clr[i][j][k] = axistmp[i][j][k];
    }
  #if 1
  } while(linferr>reconstruction_tolerance&&iter<reconstruction_maxiter);
  #else
  } while(false);
  #endif
#endif

  if(iter<reconstruction_maxiter) {
    if(verbose) {
      boil::oout<<"VOFaxisym::reconstruct_geometry converged after "
                <<iter<<" iterations! "
                <<"Final error: "<<time->current_time()<<" "<<linferr<<boil::endl;
    }
  } else {
    boil::oout<<"VOFaxisym::reconstruct_geometry did not converge after "
              <<iter<<" iterations! "
              <<"Final error: "<<time->current_time()<<" "<<linferr<<boil::endl;
    //boil::plot->plot(scp,clr,nx,nz,nalpha,"phi-clr-nx-nz-nalp", time->current_step());
    //exit(0);
  }

  return;
}
