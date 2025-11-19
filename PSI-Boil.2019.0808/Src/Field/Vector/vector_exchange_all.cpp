#include "vector.h"

/******************************************************************************/
void Vector::exchange_all(const Comp comp, const int dir) {
	
  /* through all vector components */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange_all(dir);

  /* average velocities at processor interfaces 
     should a specialized version of exchange_avg_all be written? */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange_all_avg(dir);
}
