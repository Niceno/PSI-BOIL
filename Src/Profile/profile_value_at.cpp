#include "profile.h"

/******************************************************************************/
real Profile::value_at(const real & x) const {

  assert( coord.size() == value.size() );

  real result;

  for(int i=0; i<value.size()-1; i++) {
   
    if(x >= coord[i] && x <= coord[i+1] ) {

      real w = (coord[i+1] - x) / ( coord[i+1]-coord[i] );
      
      return value[i] * w + value[i+1] * (1.0-w);
    }
  }

  /* if fails to find it, return 0 */
  return 0.0;
}

/*-----------------------------------------------------------------------------+
 '$Id: profile_value_at.cpp,v 1.2 2014/08/06 08:22:48 sato Exp $'/
+-----------------------------------------------------------------------------*/
