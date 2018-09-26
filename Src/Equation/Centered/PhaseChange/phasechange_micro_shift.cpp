#include "phasechange.h"
#include <iomanip>
#include "../../../Parallel/Out/out.h"
using namespace std;

/******************************************************************************/
void PhaseChange::micro_shift() {
/***************************************************************************//**
*  \brief advect phi in normal direction (txl, tyl, tzl)
*         tempolary: stmp
*******************************************************************************/

  for_ijk(i,j,k){
    phi[i][j][k] = dV(i,j,k) * phi[i][j][k];
  }
  phi.exchange();

#if 0
  real phisum=0.0;
  for_ijk(i,j,k){
    phisum += phi[i][j][k];
  }
  boil::cart.sum_real(&phisum);

  std::cout.setf(std::ios_base::scientific);
  std::cout<< setprecision(16);
  boil::oout<<"micro_shift: phisum-bef= "<<phisum<<"\n";
#endif
#if 0
  boil::plot->plot(phi,txl,tyl,tzl,"mdot-txl-tyl-tzl-bef",time->current_step());
#endif

#if 1
  /*-------------------------------------------+
  |  advect phi                                |
  |  time scheme: euler explicit               |
  |  advection scheme: 1st order upwind        |
  +-------------------------------------------*/
  delta=0.0;
  for(int it=0; it<6; it++){
    real dtau = 0.5*dxmin;
    //real dtau = 0.25*dxmin;
    real flux;
    for_ijk(i,j,k) {
      flux = -0.5*(txl[i][j][k]+txl[i-1][j][k])*dSx(i,j,k);
      delta[i][j][k] =-(0.5*(phi[i][j][k]+phi[i-1][j][k])*flux
                       +0.5*(phi[i][j][k]-phi[i-1][j][k])*fabs(flux));
      flux =  0.5*(txl[i][j][k]+txl[i+1][j][k])*dSx(i,j,k);
      delta[i][j][k] -= 0.5*(phi[i][j][k]+phi[i+1][j][k])*flux
                       +0.5*(phi[i][j][k]-phi[i+1][j][k])*fabs(flux);
      flux = -0.5*(tyl[i][j][k]+tyl[i][j-1][k])*dSy(i,j,k);
      delta[i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j-1][k])*flux
                       +0.5*(phi[i][j][k]-phi[i][j-1][k])*fabs(flux);
      flux =  0.5*(tyl[i][j][k]+tyl[i][j+1][k])*dSy(i,j,k);
      delta[i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j+1][k])*flux
                       +0.5*(phi[i][j][k]-phi[i][j+1][k])*fabs(flux);
      flux = -0.5*(tzl[i][j][k]+tzl[i][j][k-1])*dSz(i,j,k);
      delta[i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j][k-1])*flux
                       +0.5*(phi[i][j][k]-phi[i][j][k-1])*fabs(flux);
      flux =  0.5*(tzl[i][j][k]+tzl[i][j][k+1])*dSz(i,j,k);
      delta[i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j][k+1])*flux
                       +0.5*(phi[i][j][k]-phi[i][j][k+1])*fabs(flux);
    }
    for_ijk(i,j,k) {
      phi[i][j][k] += dtau / dV(i,j,k) * delta[i][j][k];
    }
    phi.exchange();
  }
#else
  /*-------------------------------------------+
  |  advect phi                                |
  |  time scheme: euler implicit               |
  |  advection scheme: 1st order upwind        |
  +-------------------------------------------*/
  real dtau = 1.0*dxmin;
  for(int it=0; it<1; it++){
    /* steady state part */
    for_ijk(i,j,k) {
      real flux;
      flux = -0.5*(txl[i][j][k]+txl[i-1][j][k])*dSx(i,j,k);
      stmp [i][j][k] =-(0.5*(phi[i][j][k]+phi[i-1][j][k])*flux
                       +0.5*(phi[i][j][k]-phi[i-1][j][k])*fabs(flux));
      flux =  0.5*(txl[i][j][k]+txl[i+1][j][k])*dSx(i,j,k);
      stmp [i][j][k] -= 0.5*(phi[i][j][k]+phi[i+1][j][k])*flux
                       +0.5*(phi[i][j][k]-phi[i+1][j][k])*fabs(flux);
      flux = -0.5*(tyl[i][j][k]+tyl[i][j-1][k])*dSy(i,j,k);
      stmp [i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j-1][k])*flux
                       +0.5*(phi[i][j][k]-phi[i][j-1][k])*fabs(flux);
      flux =  0.5*(tyl[i][j][k]+tyl[i][j+1][k])*dSy(i,j,k);
      stmp [i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j+1][k])*flux
                       +0.5*(phi[i][j][k]-phi[i][j+1][k])*fabs(flux);
      flux = -0.5*(tzl[i][j][k]+tzl[i][j][k-1])*dSz(i,j,k);
      stmp [i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j][k-1])*flux
                       +0.5*(phi[i][j][k]-phi[i][j][k-1])*fabs(flux);
      flux =  0.5*(tzl[i][j][k]+tzl[i][j][k+1])*dSz(i,j,k);
      stmp [i][j][k] -= 0.5*(phi[i][j][k]+phi[i][j][k+1])*flux
                       +0.5*(phi[i][j][k]-phi[i][j][k+1])*fabs(flux);
    }
    /* solve inner loop */
    delta=0.0;
    for(int loop=0; loop<8; loop++){
      /* symmetric loop */
      int ist,ied,iinc;
      int jst,jed,jinc;
      int kst,ked,kinc;
      if(loop%2==0){ist=si();ied=ei();iinc=1;}else{ist=ei();ied=si();iinc=-1;}
      if(loop%2==0){jst=sj();jed=ej();jinc=1;}else{jst=ej();jed=sj();jinc=-1;}
      if(loop%2==0){kst=sk();ked=ek();kinc=1;}else{kst=ek();ked=sk();kinc=-1;}
      std::cout<<loop<<" "<<ist<<" "<<jst<<" "<<kst<<"\n";
      for(int i=ist; i<=ied; i+=iinc){
      for(int j=jst; j<=jed; j+=jinc){
      for(int k=kst; k<=ked; k+=kinc){
      //for_ijk(i,j,k) {
        real diag = dV(i,j,k) / dtau;
        real flux,rhs;
        // west
        flux = -0.5*(txl[i][j][k]+txl[i-1][j][k])*dSx(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs  = -0.5*(flux-fabs(flux))*delta[i-1][j][k];
        // east
        flux =  0.5*(txl[i][j][k]+txl[i+1][j][k])*dSx(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs -=  0.5*(flux-fabs(flux))*delta[i+1][j][k];
        // south
        flux = -0.5*(tyl[i][j][k]+tyl[i][j-1][k])*dSy(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs -=  0.5*(flux-fabs(flux))*delta[i][j-1][k];
        // north
        flux =  0.5*(tyl[i][j][k]+tyl[i][j+1][k])*dSy(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs -=  0.5*(flux-fabs(flux))*delta[i][j+1][k];
        // bottom
        flux = -0.5*(tzl[i][j][k]+tzl[i][j][k-1])*dSz(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs -=  0.5*(flux-fabs(flux))*delta[i][j][k-1];
        // top
        flux =  0.5*(tzl[i][j][k]+tzl[i][j][k+1])*dSz(i,j,k);
        diag += 0.5*(flux+fabs(flux));
        rhs -=  0.5*(flux-fabs(flux))*delta[i][j][k+1];
        delta[i][j][k] = (stmp[i][j][k]+rhs)/diag;
      }}}
    }
    /* update */
    for_ijk(i,j,k) {
      phi[i][j][k] += delta[i][j][k];
    }
    phi.exchange();
  }

#endif

#if 0
  phisum=0.0;
  for_ijk(i,j,k){
    phisum += phi[i][j][k];
  }
  boil::cart.sum_real(&phisum);

  std::cout.setf(std::ios_base::scientific);
  std::cout<< setprecision(16);
  boil::oout<<"micro_shift: phisum-aft= "<<phisum<<"\n";
#endif

#if 0
  boil::plot->plot(phi,txl,tyl,tzl,"mdot-txl-tyl-tzl-aft",time->current_step());
  exit(0);
#endif

  for_ijk(i,j,k){
    phi[i][j][k] = phi[i][j][k] / dV(i,j,k);
  }
  phi.exchange();

  //mdot_cut();
  //phi.exchange();

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: phasechange_micro_shift.cpp,v 1.1 2014/08/06 08:19:37 sato Exp $'/
+-----------------------------------------------------------------------------*/
