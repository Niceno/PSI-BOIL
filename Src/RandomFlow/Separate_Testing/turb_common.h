#include <iostream>
#include <sstream>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>

#define DIM 3

#include "../Global/global_malloc.h"

real sclp(const real * A, const real * B); 
void vecp(real * A, const real * B, const real * C); 
void seed(); 
real rnd(); 
real gauss();
void gaussn(real **Y, real d, int n, int m); 
void gaussn(real *Y, real d, int n); 

/*-----------------------------------------------------------------------------+
 '$Id: turb_common.h,v 1.2 2008/10/21 11:54:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
