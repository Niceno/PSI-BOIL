#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::dirac_source_terms() {
/***************************************************************************//**
*  \brief Corrects heat flux at solid-fluid boundaries for extrapolation.
******************************************************************************/

  /*--------------+
  | immersed body |
  +--------------*/
  if(solid()) {
     for_ijk(i,j,k) {
       if(dom->ibody().off(i,j,k)) {
         if(dom->ibody().on(i-1,j,k)) {
           txv[i][j][k] += cht.dirac_wall_source(i-1,j,k);
           txl[i][j][k] += cht.dirac_wall_source(i-1,j,k);
         }
         if(dom->ibody().on(i+1,j,k)) {
           txv[i][j][k] -= cht.dirac_wall_source(i+1,j,k);
           txl[i][j][k] -= cht.dirac_wall_source(i+1,j,k);
         }
         if(dom->ibody().on(i,j-1,k)) {
           tyv[i][j][k] += cht.dirac_wall_source(i,j-1,k);
           tyl[i][j][k] += cht.dirac_wall_source(i,j-1,k);
         }
         if(dom->ibody().on(i,j+1,k)) {
           tyv[i][j][k] -= cht.dirac_wall_source(i,j+1,k);
           tyl[i][j][k] -= cht.dirac_wall_source(i,j+1,k);
         }
         if(dom->ibody().on(i,j,k-1)) {
           tzv[i][j][k] += cht.dirac_wall_source(i,j,k-1);
           tzl[i][j][k] += cht.dirac_wall_source(i,j,k-1);
         }
         if(dom->ibody().on(i,j,k+1)) {
           tzv[i][j][k] -= cht.dirac_wall_source(i,j,k+1);
           tzl[i][j][k] -= cht.dirac_wall_source(i,j,k+1);
         }
       } /* center off */
     } /* ijk */
  } /* solid exists */

  return;
}
