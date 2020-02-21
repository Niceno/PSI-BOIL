#include "vof.h"

/******************************************************************************/
void VOF::advance_x(Scalar & scp) {
  
  // advance in the x-direction

  Comp m = Comp::u();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei()+1; i++)
  //for(int j = sj(); j <=ej(); j++)
  //for(int k = sk(); k <=ek(); k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
    Sign sig = Sign::neg();

    // upwind i-index
    int iup = i-1;
    if((*u)[m][i][j][k]<0.0) iup = i;             

    real dt = time->dt();

    if (scp[iup][j][k] < boil::pico) {

      f = scp[iup][j][k] * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(scp[iup][j][k]>1.0-boil::pico) {

      f = scp[iup][j][k] * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dxc(iup)==0.0 ) {

        f = scp[iup][j][k] * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/scp.dxc(iup);
        real scpup = scp[iup][j][k];
        f = dV(iup,j,k)*calc_flux(g,scpup,nx[iup][j][k],
                                          ny[iup][j][k],
                                          nz[iup][j][k],
                                          nalpha[iup][j][k]);
      }
    }

    // update stmp
    stmp[i-1][j][k] = stmp[i-1][j][k] - f;
    stmp[i  ][j][k] = stmp[i  ][j][k] + f;
    vflow[m][i][j][k] = f;

  }

  return;

}

