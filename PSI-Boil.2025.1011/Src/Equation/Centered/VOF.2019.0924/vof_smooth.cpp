#include "vof.h"
#include <iomanip>

/******************************************************************************/
void VOF::smooth(const Scalar & sca, Scalar & scb, const int itnum) {
/***************************************************************************//**
* \brief  Smooth/sharpen color function and cut-off.
*         If itnum=0, then only the cut-off function works.
*            input    : sca, itnum
*            output   : scb
*            temporary: dflag
*******************************************************************************/

  for_aijk(i,j,k)
    scb[i][j][k]=sca[i][j][k];

  if (itnum>=1) {
#if 0
    /*----------+
    |   Yabe    |
    +----------*/
    /* coefficients for smooth */
    const real c1 = 1.0/(6.0+12.0/sqrt(2.0)+8/sqrt(3.0));
    const real c2 = c1/sqrt(2.0);
    const real c3 = c1/sqrt(3.0);

    /* iterate */
    for(int it=0; it<itnum; it++) {
      for_ijk(i,j,k) {
        stmp[i][j][k] =0.5*scb[i][j][k] + 0.5/(1.0+6.0*c1+12.0*c2+8.0*c3)
                         *(scb[i][j][k]
                      +c1*(scb[i-1][j][k]+scb[i+1][j][k]
                          +scb[i][j-1][k]+scb[i][j+1][k]
                          +scb[i][j][k-1]+scb[i][j][k+1])
                      +c2*(scb[i-1][j-1][k]+scb[i-1][j+1][k] 
                          +scb[i+1][j-1][k]+scb[i+1][j+1][k]
                          +scb[i-1][j][k-1]+scb[i-1][j][k+1]
                          +scb[i+1][j][k-1]+scb[i+1][j][k+1]
                          +scb[i][j-1][k-1]+scb[i][j-1][k+1]
                          +scb[i][j+1][k-1]+scb[i][j+1][k+1])
                      +c3*(scb[i-1][j-1][k-1]+scb[i-1][j-1][k+1]
                          +scb[i-1][j+1][k-1]+scb[i-1][j+1][k+1]
                          +scb[i+1][j-1][k-1]+scb[i+1][j-1][k+1]
                          +scb[i+1][j+1][k-1]+scb[i+1][j+1][k+1]));
      }
      //insert_bc(dflag);
      //dflag.exchange_all();

      for_ijk(i,j,k)
        scb[i][j][k]=stmp[i][j][k];

      scb.bnd_update(); // boundary condition for scb
      scb.exchange();
    }
#endif
#if 1
    /*-----------------------------+
    |  diffusion, explicit Jakobi  |
    +-----------------------------*/
    real dtau=0.125;
    /* iterate */
    for(int it=0; it<itnum; it++) {
      for_ijk(i,j,k) {
        stmp[i][j][k] = scb[i-1][j][k] -2.0*scb[i][j][k] +scb[i+1][j][k]
                      + scb[i][j-1][k] -2.0*scb[i][j][k] +scb[i][j+1][k]
                      + scb[i][j][k-1] -2.0*scb[i][j][k] +scb[i][j][k+1];
      }
      //insert_bc(nmag);
      //nmag.exchange_all();

      // update
      //for_aijk(i,j,k){
      for_ijk(i,j,k){
        //real coef=std::min(1.0,pow(abs(2.0*sca[i][j][k]-1.0),1.0));
        //scb[i][j][k]=scb[i][j][k]+dtau*coef*nmag[i][j][k];
        scb[i][j][k]=scb[i][j][k]+dtau*stmp[i][j][k];
      }
      //insert_bc_diffeq(scb); // boundary condition for scb
      scb.bnd_update(); // boundary condition for scb
      scb.exchange();
    }

#endif
  }

  /*----------+
  |  cut-off  |
  +----------*/
  for_aijk(i,j,k)
    scb[i][j][k]=std::max(0.0,(std::min(1.0,scb[i][j][k])));

#if 0
  boil::plot->plot(sca,scb, "sca-scb", time->current_step());
  exit(0);
#endif

  return;
}
