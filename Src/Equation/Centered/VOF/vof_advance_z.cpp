#include "vof.h"

/******************************************************************************/
void VOF::advance_z(Scalar & scp) {

  /* advance in the z-direction */

  Comp m = Comp::w();

  //for_vmijk((*u),m,i,j,k){
  //for(int i = si(); i <=ei(); i++)
  //for(int j = sj(); j <=ej(); j++)
  //for(int k = sk(); k <=ek()+1; k++) {
  for_wvmijk((*u),m,i,j,k) { /* this expands to the three lines above */

    /* flux */
    real f;

    // upwind j-index
    int kup = k-1;
    if((*u)[m][i][j][k]<0.0) kup = k;             

    real dt=time->dt();

    if (scp[i][j][kup] < boil::pico) {

      f = scp[i][j][kup] * dSz(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(scp[i][j][kup]>1.0-boil::pico) {

      f = scp[i][j][kup] * dSz(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (scp.dzc(kup)==0.0) {

        f = scp[i][j][kup] * dSz(i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/scp.dzc(kup);
#if 1
        real scpup = scp[i][j][kup];
        f = dV(i,j,kup)*calc_flux(g,scpup,nz[i][j][kup],
                                          ny[i][j][kup],
                                          nx[i][j][kup]);
#else
        if (g==0.0) {
          f = 0.0;
        } else {
          // color function upwind
          real c = scp[i][j][kup];
      
          // calculate vn1, vn2, vn3: normal vector at face center
          real vn1 = -nx[i][j][kup];
          real vn2 = -ny[i][j][kup];
          real vn3 = -nz[i][j][kup];

          real absg = fabs(g);
          real vm1 = fabs(vn1);
          real vm2 = fabs(vn2);
          real vm3 = fabs(vn3)+boil::pico;
          real qa = 1.0/(vm1+vm2+vm3);
          vm1 *= qa;
          vm2 *= qa;
          vm3 *= qa;
          real alpha = calc_alpha(c, vm1, vm2, vm3);
              
          real ra = vm3 * (1.0 - absg);
          qa = 1.0/(1.0-ra);
          if (g*vn3 > 0) alpha = alpha -ra;
          vm3 = vm3 * absg;
      
          // calculate f: flux
          f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g * dV(i,j,kup);
        }
#endif
      }
    }

    // update stmp
    stmp[i][j][k-1] = stmp[i][j][k-1] - f;
    stmp[i  ][j][k] = stmp[i  ][j][k] + f;

#if 0
    if((k==100||k==101)&&j==3&&i==100) {
      std::cout<<"advance_z:"<<f<<" "<<k<<" "<<scp[i][j][k-1]<<"\n";
    }
#endif

  }

  return;

}

