#include "vof.h"

/******************************************************************************/
void VOF::flood(Scalar & scp,const real mult) {
/***************************************************************************//**
*  \brief Flood walls with either 0 (mult<0) or 1 (mult>0).
*******************************************************************************/

  real colref = mult < 0. ? 0. : 1.;

  for_avijk(scp,i,j,k) {
    if(   dom->ibody().off(i,j,k) 
       || (i<si() && iminw) || (i>ei() && imaxw)
       || (j<sj() && jminw) || (j>ej() && jmaxw)
       || (k<sk() && kminw) || (k>ek() && kmaxw)
      ) {
      scp[i][j][k] = colref;
    }
  }

  return;
}

