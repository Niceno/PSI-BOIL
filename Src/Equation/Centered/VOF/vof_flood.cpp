#include "vof.h"

/******************************************************************************/
void VOF::flood(Scalar & scp,const real mult) {
/***************************************************************************//**
*  \brief Flood walls with either 0 (mult<0) or 1 (mult>0).
*******************************************************************************/

  real colref = mult < 0. ? 0. : 1.;

  for_avijk(scp,i,j,k) {
    if(   dom->ibody().off(i,j,k) 
       || (i<si() && bflag_struct.iminw) || (i>ei() && bflag_struct.imaxw)
       || (j<sj() && bflag_struct.jminw) || (j>ej() && bflag_struct.jmaxw)
       || (k<sk() && bflag_struct.kminw) || (k>ek() && bflag_struct.kmaxw)
      ) {
      scp[i][j][k] = colref;
    }
  }

  return;
}

