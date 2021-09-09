#include "vof.h"

/******************************************************************************/
void VOF::advance_z(const Scalar & scp, Scalar & cellvol) {

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

    // upwind k-index
    int kup = k-1;
    if((*u)[m][i][j][k]<0.0) kup = k;             

    real dt=time->dt();

    if (scp[i][j][kup] < boil::pico) {

      f = scp[i][j][kup] * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(scp[i][j][kup]>1.0-boil::pico) {

      f = scp[i][j][kup] * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dzc(kup)==0.0) {

        f = scp[i][j][kup] * dSz(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/scp.dzc(kup);
        real scpup = scp[i][j][kup];
        f = dV(i,j,kup)*calc_flux(g,scpup,nz[i][j][kup],
                                          ny[i][j][kup],
                                          nx[i][j][kup],
                                          nalpha[i][j][kup]);
      }
    }

    /* update cellvol */
    cellvol[i][j][k-1] = cellvol[i][j][k-1] - f;
    cellvol[i  ][j][k] = cellvol[i  ][j][k] + f;
    vflow[m][i][j][k] = f;

  }

  return;

}

