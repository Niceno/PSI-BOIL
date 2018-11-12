#include "vector.h"

/******************************************************************************/
void Vector::exchange(const Comp comp, const int dir) {
	
  /* through all vector components */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange(dir);

  /* average velocities at processor interfaces */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange_avg(dir);

}

/******************************************************************************/
void Vector::exchange(const int * i, const Comp comp, const int dir) {

  /* through all vector components */
  for_m(m)
    if(comp == Comp::undefined() || comp == m)
      vec[m].exchange(i, dir);

}
