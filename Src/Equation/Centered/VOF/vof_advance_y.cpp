#include "vof.h"

/******************************************************************************/
void VOF::advance_y(const Scalar & scp) {

  /* advance in the y-direction */

  Comp m = Comp::v();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej()+1; j++)
  //for(int k = sk(); k <=ek(); k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;
    Sign sig = Sign::neg();

    // upwind j-index
    int jup = j-1;
    if((*u)[m][i][j][k]<0.0) jup = j;             

    real dt=time->dt();

    if (scp[i][jup][k]<boil::pico) {

      f = scp[i][jup][k] * dSy(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(scp[i][jup][k]>1.0-boil::pico) {

      f = scp[i][jup][k] * dSy(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dyc(jup)==0.0) {

        f = scp[i][jup][k] * dSy(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/scp.dyc(jup);
        real scpup = scp[i][jup][k];
        f = dV(i,jup,k)*calc_flux(g,scpup,ny[i][jup][k],
                                          nx[i][jup][k],
                                          nz[i][jup][k],
                                          nalpha[i][jup][k]);
      }
    }

    // update stmp
    stmp[i][j-1][k] = stmp[i][j-1][k] - f;
    stmp[i][j  ][k] = stmp[i][j  ][k] + f;
    vflow[m][i][j][k] = f;

  }

  return;

}

