#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::cal_gradt() {
/***************************************************************************//*** 
*  \brief calculate temperature gradient in all cells  
*******************************************************************************/

#if 0
  *  iflag = -3 : color function <  clrsurf & no interface nearby
*  iflag = +3 : color function >= clrsurf & no interface nearby
*  iflag = -2 : color function <  clrsurf & interface nearby
*  iflag = +2 : color function >= clrsurf & interface nearby
*  iflag = -1 : color function <  clrsurf & at interface
*  iflag = +1 : color function >= clrsurf & at interface
*  iflag = -1000 solid cell
#endif

  for_ijk(i,j,k) {
    const bool is_solid = dom->ibody().off(i,j,k);
    real gt = gradt1D(is_solid,Comp::i(),i,j,k);
  }
}
