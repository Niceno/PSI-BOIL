#include "vof.h"

/******************************************************************************/
void VOF::advance_y(Scalar & scp) {

  /* advance in the y-direction */

  Comp m = Comp::v();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej()+1; j++)
  //for(int k = sk(); k <=ek(); k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;

    // upwind j-index
    int jup = j-1;
    if((*u)[m][i][j][k]<0.0) jup = j;             

    real dt=time->dt();

    if (scp[i][jup][k]<boil::pico) {

      f = scp[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(scp[i][jup][k]>1.0-boil::pico) {

      f = scp[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dyc(jup)==0.0) {

        f = scp[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/scp.dyc(jup);
#if 1
        real scpup = scp[i][jup][k];
        f = dV(i,jup,k)*calc_flux(g,scpup,ny[i][jup][k],
                                          nx[i][jup][k],
                                          nz[i][jup][k]);
#else 
        if (g==0.0) {
          f = 0.0;
        } else {
          // color function upwind
          real c = scp[i][jup][k];

          // calculate vn1, vn2, vn3: normal vector at face center
          real vn1 = -nx[i][jup][k];
          real vn2 = -ny[i][jup][k];
          real vn3 = -nz[i][jup][k];

          real absg = fabs(g);
          real vm1 = fabs(vn1);
          real vm2 = fabs(vn2);
          real vm3 = fabs(vn3)+boil::pico;
          real qa = 1.0/(vm1+vm2+vm3);
          vm1 *= qa;
          vm2 *= qa;
          vm3 *= qa;
          real alpha = calc_alpha(c, vm1, vm2, vm3);
      
          real ra = vm2 * (1.0 - absg);
          qa = 1.0/(1.0-ra);
          if (g*vn2 > 0) alpha = alpha -ra;
          vm2 = vm2 * absg;

          // calculate f: flux
          f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g * dV(i,jup,k);
        }
#endif
      }
    }

    // update stmp
    stmp[i][j-1][k] = stmp[i][j-1][k] - f;
    stmp[i][j  ][k] = stmp[i][j  ][k] + f;
#if 0
    if(i==100&&j==3&&k==100) {
      std::cout<<"advance_y:"<<f<<" "<<j<<" "<<scp[i][j][k]<<"\n";
    }
#endif

  }

  return;

}

