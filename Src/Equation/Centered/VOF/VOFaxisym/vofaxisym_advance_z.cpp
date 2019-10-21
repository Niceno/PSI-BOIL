#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::advance_z(Scalar & scp) {

  /* advance in the z-direction */

  Comp m = Comp::w();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej(); j++)
  //for(int k = sk(); k <=ek()+1; k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
    Sign sig = Sign::neg();

    /* upwind k-index */
    int kup = k-1;
    if((*u)[m][i][j][k]<0.0) kup = k;             

    real dt=time->dt();
#if 1
    real refval = clr[i][j][kup];
#else
    real refval = phi[i][j][kup];
#endif

    if (refval < boil::pico) {

      f = refval * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(refval>1.0-boil::pico) {

      f = refval * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dzc(kup)==0.0) {

        f = refval * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else { 

        /* calculate g: CFL upwind */
        real g = ((*u)[m][i][j][k])*dt/scp.dzc(kup);
#if 1
        f = dV(i,j,kup)*calc_flux_axisymmetric(g,i,j,kup,Comp::w());
#else
        real scpup = scp[i][j][kup];
        f = dV(i,j,kup)*calc_flux(g,scpup,nz[i][j][kup],
                                          ny[i][j][kup],
                                          nx[i][j][kup]);
#endif
      }
    }

    /* update stmp */
    stmp[i][j][k-1] = stmp[i][j][k-1] - f;
    stmp[i  ][j][k] = stmp[i  ][j][k] + f;

  }

  return;

}

