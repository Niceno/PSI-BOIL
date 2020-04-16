#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::cal_hf(const Scalar * diff_eddy) {
/***************************************************************************//*** 
*  \brief calculate heat flux in all cells  
*******************************************************************************/

  for_ijk(i,j,k) {
    const bool is_solid = dom->ibody().off(i,j,k);
    const real lmb = lambda(i,j,k,diff_eddy);

    real gtx = gradt1D(is_solid,Comp::i(),i,j,k);
    real gty = gradt1D(is_solid,Comp::j(),i,j,k);
    real gtz = gradt1D(is_solid,Comp::k(),i,j,k);

    /* ghost values */
    if(is_solid) {
      txv[i][j][k] = txl[i][j][k] = lmb*gtx;
      tyv[i][j][k] = tyl[i][j][k] = lmb*gty;
      tzv[i][j][k] = tzl[i][j][k] = lmb*gtz;
    } else {
      /* vapor */
      if(clr[i][j][k]<clrsurf) {
        txv[i][j][k] = lmb*gtx;
        tyv[i][j][k] = lmb*gty;
        tzv[i][j][k] = lmb*gtz;
        txl[i][j][k] = 0.0;
        tyl[i][j][k] = 0.0;
        tzl[i][j][k] = 0.0;
      /* liquid */
      } else {
        txv[i][j][k] = 0.0;
        tyv[i][j][k] = 0.0;
        tzv[i][j][k] = 0.0;
        txl[i][j][k] = lmb*gtx;
        tyl[i][j][k] = lmb*gty;
        tzl[i][j][k] = lmb*gtz;
      }
    }
  } /* ijk */

  return;
}
