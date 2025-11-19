#include "vectorbool.h"

/******************************************************************************/
void VectorBool::exchange_all(const Comp comp, const int dir) {
	
  /* through all vector components */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange_all(dir);
}
