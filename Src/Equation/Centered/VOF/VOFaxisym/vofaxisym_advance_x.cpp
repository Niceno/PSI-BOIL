#include "vofaxisym.h"

/******************************************************************************/
void VOFaxisym::advance_x(const Scalar & scp, Scalar & cellvol) {
  
  /* advance in the x-direction */

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
    real refval = phi[iup][j][k];
    real n2sum = nx[iup][j][k]*nx[iup][j][k]+nz[iup][j][k]*nz[iup][j][k];

    if (refval < boil::pico || n2sum<0.5) {

      f = refval * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(refval>1.0-boil::pico) {

      f = refval * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else { 

      if (scp.dxc(iup)==0.0 ) {

        f = refval * dSx(sig,i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else { 

        /* calculate g: CFL upwind */
        real g = ((*u)[m][i][j][k])*dt/scp.dxc(iup);
#if 1
        f = dV(iup,j,k)*calc_flux_axisymmetric(g,iup,j,k,Comp::u());
#else
        real scpup = scp[iup][j][k];
        f = dV(iup,j,k)*calc_flux(g,scpup,nx[iup][j][k],
                                          ny[iup][j][k],
                                          nz[iup][j][k],
                                          nalpha[iup][j][k]);
#endif

      }
    }

    /* update cellvol */
    cellvol[i-1][j][k] = cellvol[i-1][j][k] - f;
    cellvol[i  ][j][k] = cellvol[i  ][j][k] + f;
    vflow[m][i][j][k] = f;

  }

  return;

}

