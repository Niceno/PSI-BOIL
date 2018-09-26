#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>
#include <ctime>
#include <cstdlib>

#include "global_checking.h"
#include "global_precision.h"
#include "../Ravioli/range.h"

namespace boil {

extern void random_seed();
extern void random_seed(int i);
extern real random_number(); 
extern real random_number(const Range<real> &); 

} /* boil */

#endif

/*-----------------------------------------------------------------------------+
 '$Id: global_random.h,v 1.4 2015/11/12 07:47:16 sato Exp $'/
+-----------------------------------------------------------------------------*/
