#include "vofaxisym.h"

void VOFaxisym::reconstruct_geometry() {
  reconstruct_geometry(phi);
  return;
}

/******************************************************************************/
void VOFaxisym::reconstruct_geometry(const Scalar & scp) {
/***************************************************************************//**
 \brief Reconstruct the geometry of the interface, i.e. calculate the normal
        vector, line constant alpha and color function. The reconstruction pro-
        cess is iterative.
    line: vm1*x + vm2*z = alpha
    output: nx, ny = 0.0, nz, nalpha, c
*******************************************************************************/

  real linferr(boil::unreal);
  int iter(0);
#if 0
  /* simplified iteration */

  /* set K0 = 1.0 */
  for_avijk(Ktmp,i,j,k) {
    Ktmp[i][j][k] = 1.0;
  }

  /* iteration procedure */
  do {
    /* set c = phi/K */
    for_vijk(clr,i,j,k) {
      clr[i][j][k] = scp[i][j][k]/Ktmp[i][j][k];
    }
    clr.bnd_update();
    clr.exchange_all();

    /* Ktmp is internally overwritten in this function! */
    color_to_vf(clr,axistmp);

    /* calculate Linf error norm of reconstruction in phi-space */
    linferr = linf_scalar_error(axistmp,scp);

    iter++;
    boil::oout<<"VOFaxisym::reconstruct_geometry: "<<time->current_time()
              <<" "<<iter<<" "<<linferr<<boil::endl;
  #if 1
  } while(linferr>reconstruction_tolerance&&iter<reconstruction_maxiter);
  #else
  } while(false);
  #endif
#else
  /* direct iteration */

  /* set c0 = phi */
  for_avijk(clr,i,j,k) {
    clr[i][j][k] = std::max(0.0,std::min(1.0,scp[i][j][k]));
  }

  #if 1
  /* iteration procedure */
  do {
    norm_axisymmetric(clr);
    backward_axisymmetric(phi,nalpha);
    forward(axistmp);
    #if 0
    /* calculate Linf error norm of reconstruction in c-space */
    linferr = linf_scalar_error(axistmp,clr);
    #else
    color_to_vf(axistmp,axistmp2);
    /* calculate Linf error norm of reconstruction in phi-space */
    linferr = linf_scalar_error(axistmp2,phi);
    #endif

    iter++;
    boil::oout<<"VOFaxisym::reconstruct_geometry: "<<time->current_time()
              <<" "<<iter<<" "<<linferr<<boil::endl;

    for_avijk(clr,i,j,k) {
      clr[i][j][k] = axistmp[i][j][k];
    }
    #if 1
  } while(linferr>reconstruction_tolerance&&iter<reconstruction_maxiter);
    #else
  } while(false);
    #endif
  #else
  norm_axisymmetric(clr);
  extract_alpha(clr);
  #endif
#endif

  if(iter<reconstruction_maxiter) {
    boil::oout<<"VOFaxisym::reconstruct_geometry converged after "
              <<iter<<" iterations! "
              <<"Final error: "<<linferr<<boil::endl;
  } else {
    boil::oout<<"VOFaxisym::reconstruct_geometry did not converge after "
              <<iter<<" iterations! "
              <<"Final error: "<<linferr<<boil::endl;
    //boil::plot->plot(phi,clr,nx,nz,nalpha,"phi-clr-nx-nz-nalp", time->current_step());
    //exit(0);
  }

  return;
}
