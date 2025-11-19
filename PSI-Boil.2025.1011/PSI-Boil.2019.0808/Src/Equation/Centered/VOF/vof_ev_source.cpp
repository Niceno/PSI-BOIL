#include "vof.h"

//#define DEBUG

/******************************************************************************/
void VOF::ev_calculate_source(const ScalarInt & pflag, const Vector & u,
                              Scalar & fext) {
/***************************************************************************//**
*  \brief Calculate source term from divergence of velocity on sol. domain.
*  pflag = -2,-1,1,-2 : solution domain
*  fext = div u / dt
*******************************************************************************/

  /* initialize fext */
  fext = 0.;

  /* source calculation */
  for_ijk(i,j,k) {
    if(abs(pflag[i][j][k])<3) {
      fext[i][j][k] = u.outflow(i,j,k)*time->dti();
    }
  }
  fext.exchange();

  return;
}
