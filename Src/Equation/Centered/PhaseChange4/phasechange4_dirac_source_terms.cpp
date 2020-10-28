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
           txv[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
           txl[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
         }
         if(dom->ibody().on(i+1,j,k)) {
           txv[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
           txl[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
         }
         if(dom->ibody().on(i,j-1,k)) {
           tyv[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
           tyl[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
         }
         if(dom->ibody().on(i,j+1,k)) {
           tyv[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
           tyl[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
         }
         if(dom->ibody().on(i,j,k-1)) {
           tzv[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
           tzl[i][j][k] += cht.heat_transfer_wall_model().dirac_wall_source;
         }
         if(dom->ibody().on(i,j,k+1)) {
           tzv[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
           tzl[i][j][k] -= cht.heat_transfer_wall_model().dirac_wall_source;
         }
       } /* center off */
     } /* ijk */
  } /* solid exists */

  return;
}
