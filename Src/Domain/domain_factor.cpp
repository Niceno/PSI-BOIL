#include "domain.h"

/******************************************************************************/
void Domain::factor(int n, int * factor, int * number) const {
 
  (*number) = 0;
  
  int m=n-1;
  for(;;) {
    if(n % m == 0) {
      factor[(*number)++] = n/m;
      n = m;
    }
    m--;	  
    if(m == 0) break;
    if(m == 1) {
      factor[(*number)++] = n;
      break;
    }
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: domain_factor.cpp,v 1.4 2009/03/14 11:37:58 niceno Exp $'/
+-----------------------------------------------------------------------------*/
