#include "global_random.h"

namespace boil {

/******************************************************************************/
void random_seed() {
  const int sequences=10000000;
  srandom(time(0) % sequences);
}

/******************************************************************************/
void random_seed(int i) {
  srandom(i);
}

/******************************************************************************/
real random_number() {
/*----------------------------------------------+
|  Returns a random number between 0.0 and 1.0  |
+----------------------------------------------*/
  return (real)random()/RAND_MAX;
}

/******************************************************************************/
real random_number(const Range<real> & rng) {
/*-------------------------------------------------+
|  Returns a random number in the specified range  |
+-------------------------------------------------*/

  real length = rng.last() - rng.first();
  
  assert(length > 0);
  
  return rng.first() + random_number() * length;
}

} /* boil */
