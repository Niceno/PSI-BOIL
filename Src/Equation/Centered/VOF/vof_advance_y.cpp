#include "vof.h"

/******************************************************************************/
void VOF::advance_y() {

  // advance in the y-direction

  Comp m = Comp::v();

  for(int i = si(); i <=ei(); i++)
  for(int j = sk(); j <=ej()+1; j++)
  for(int k = sj(); k <=ek(); k++) {

    // flux
    real f;

    // upwind j-index
    int jup = j-1;
    if((*u)[m][i][j][k]<0.0) jup = j;             

    real dt=time->dt();

    if (phi[i][jup][k]<boil::pico) {

      f = phi[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(phi[i][jup][k]>1.0-boil::pico) {

      f = phi[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (phi.dyc(jup)==0.0) {

        f = phi[i][jup][k] * dSy(i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/phi.dyc(jup);

        if (g==0.0) {
          f = 0.0;
        } else {
          // color function upwind
          real c = phi[i][jup][k];

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
      }
    }

    // update stmp
    stmp[i][j-1][k] = stmp[i][j-1][k] - f;
    stmp[i][j  ][k] = stmp[i][j  ][k] + f;
#if 0
    if(i==100&&j==3&&k==100) {
      std::cout<<"advance_y:"<<f<<" "<<j<<" "<<phi[i][j][k]<<"\n";
    }
#endif

  }

  return;

}

