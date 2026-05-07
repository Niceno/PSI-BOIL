#include "vof.h"
#include <iomanip>
using namespace boil;

/******************************************************************************/
void VOF::smooth(const Scalar & sca, Scalar & scb, const int itnum) {
/***************************************************************************//**
* \brief  Smooth/sharpen color function and cut-off.
*         If itnum=0, then only the cut-off function works.
*            input    : sca, itnum
*            output   : scb
*            temporary: dflag
*******************************************************************************/

  // set initial value
  for_aijk(i,j,k)
    scb[i][j][k]=sca[i][j][k];

  /*-----------------------------+
  |  diffusion, explicit Jakobi  |
  +-----------------------------*/
  /* iterate */
  for(int it=0; it<itnum; it++) {
    real diff=0.0;
    for_ijk(i,j,k) {
      real lambda = 1.0*std::max(phi.dxc(i),std::max(phi.dyc(j),phi.dzc(k)));
      real dtau=1.0;
      real dxm = phi.dxw(i);
      real dxp = phi.dxe(i);
      real dym = phi.dys(j);
      real dyp = phi.dyn(j);
      real dzm = phi.dzb(k);
      real dzp = phi.dzt(k);
      real diag = 1.0 + dtau * lambda *
                ( 2.0/(dxm+dxp)*(1.0/dxp+1.0/dxm)
                + 2.0/(dym+dyp)*(1.0/dyp+1.0/dym)
                + 2.0/(dzm+dzp)*(1.0/dzp+1.0/dzm));
      real rhs = scb[i][j][k] + dtau * lambda *
                ( 2.0/(dxm+dxp)*(scb[i+1][j][k]/dxp+scb[i-1][j][k]/dxm)
                + 2.0/(dym+dyp)*(scb[i][j+1][k]/dyp+scb[i][j-1][k]/dym)
                + 2.0/(dzm+dzp)*(scb[i][j][k+1]/dzp+scb[i][j][k-1]/dzm));
      stmp[i][j][k] = rhs/diag;
      diff += pow(scb[i][j][k]-stmp[i][j][k],2);
      if(pow(scb[i][j][k]-stmp[i][j][k],2)>1) {
        std::cout<<"large diff:"<<i<<" "<<j<<" "<<k<<" "<<scb[i][j][k]<<" "<<stmp[i][j][k]<<"\n";
        exit(0);
      }
    }
    boil::oout<<"vof_smooth: "<<it<<" "<<sqrt(diff)<<"\n";
    stmp.bnd_update();   // BUG IN bnd_update: the range is 5<i,j,k<...! vof.cpp
    stmp.exchange_all();

    // update
    for_aijk(i,j,k){
      scb[i][j][k]=stmp[i][j][k];
    }
  }

  /*----------+
  |  cut-off  |
  +----------*/
#if 1
  for_aijk(i,j,k)
    scb[i][j][k]=maxr(0.0,(minr(1.0,scb[i][j][k])));
#endif

#if 0
  boil::plot->plot(sca,scb, "sca-scb", time->current_step());
  exit(0);
#endif

  return;
}
