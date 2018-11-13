#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::cal_gradclr() {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr
*******************************************************************************/

  /* cell centered base, second order */
  for_ijk(i,j,k) {
      real gradient = 0.0;

      real dxm = -clr.dxw(i);
      real dxp = clr.dxe(i);
      real dym = -clr.dys(j);
      real dyp = clr.dyn(j);
      real dzm = -clr.dzb(k);
      real dzp = clr.dzt(k);
 
      real c00 = clr[i][j][k];
      real cxm = clr[i-1][j][k];
      real cxp = clr[i+1][j][k];
      real cym = clr[i][j-1][k];
      real cyp = clr[i][j+1][k];
      real czm = clr[i][j][k-1];
      real czp = clr[i][j][k+1];
 
#if 1
      real gradx = grad_2nd(c00,cxm,cxp,dxm,dxp); 
      real grady = grad_2nd(c00,cym,cyp,dym,dyp); 
      real gradz = grad_2nd(c00,czm,czp,dzm,dzp); 
#else
      real gradx = (cxp-cxm)/(dxp-dxm); 
      real grady = (cyp-cym)/(dyp-dym); 
      real gradz = (czp-czm)/(dzp-dzm); 
#endif
      gradclr[i][j][k] = sqrt(gradx*gradx+grady*grady+gradz*gradz);
  }

  real sum = 0.0;
  for_ijk(i,j,k)
    sum += gradclr[i][j][k]*dV(i,j,k);
  boil::oout<<"PCV::gradclr "<<time->current_time()<<" "<<sum<<boil::endl;

  return;
}

void PhaseChangeVOF::set_gradclr_ext(const Scalar & heavi) {

 for_ijk(i,j,k) {
   if(heavi[i][j][k]>boil::atto&&heavi[i][j][k]<1.-boil::atto)
     gradclr[i][j][k] = 1./phi.dxc(i);
   else 
     gradclr[i][j][k] = 0.;
 }

}
