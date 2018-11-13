#include "vof.h"

/******************************************************************************/
void VOF::advance_x() {

  
  // advance in the x-direction

  Comp m = Comp::u();

  //for_vmijk((*u),m,i,j,k){
  for(int i = si(); i <=ei()+1; i++)
  for(int j = sk(); j <=ej(); j++)
  for(int k = sj(); k <=ek(); k++) {

    // flux
    real f;

    // upwind i-index
    int iup = i-1;
    if((*u)[m][i][j][k]<0.0) iup = i;             

    real dt = time->dt();

    if (phi[iup][j][k] < boil::pico) {

      f = phi[iup][j][k] * dSx(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else if(phi[iup][j][k]>1.0-boil::pico) {

      f = phi[iup][j][k] * dSx(i,j,k) * ((*u)[m][i][j][k]) * dt;

    } else {

      if (phi.dxc(iup)==0.0 ) {

        f = phi[iup][j][k] * dSx(i,j,k) * ((*u)[m][i][j][k]) * dt;

      } else {

        // calculate g: CFL upwind
        real g = ((*u)[m][i][j][k])*dt/phi.dxc(iup);

        if (g==0.0) {
          f=0.0;
        } else {
          // color function upwind
          real c = phi[iup][j][k];
  
          // calculate vn1, vn2, vn3: normal vector at face center
          real vn1 = -nx[iup][j][k];
          real vn2 = -ny[iup][j][k];
          real vn3 = -nz[iup][j][k];

          real absg = fabs(g);
          real vm1 = fabs(vn1);
          real vm2 = fabs(vn2);
          real vm3 = fabs(vn3)+boil::pico;
          real qa = 1.0/(vm1+vm2+vm3);
          vm1 *= qa;
          vm2 *= qa;
          vm3 *= qa;
          real alpha = calc_alpha(c, vm1, vm2, vm3);
      
          real ra = vm1 * (1.0 - absg);
          qa = 1.0/(1.0-ra);
          if (g*vn1 > 0) alpha = alpha -ra;
          vm1 = vm1 * absg;

          // calculate f: flux
          f = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa) * g * dV(iup,j,k);
        }
      }
    }

    // update stmp
    stmp[i-1][j][k] = stmp[i-1][j][k] - f;
    stmp[i  ][j][k] = stmp[i  ][j][k] + f;

#if 0
    if((i==100||i==101)&&j==3&&k==100) {
      std::cout<<"advance_x:"<<f<<" "<<i<<" "<<phi[i-1][j][k]<<"\n";
      //if (phi[i-1][j][k]>0.1) exit(0);
    }
#endif

  }

  return;

}

