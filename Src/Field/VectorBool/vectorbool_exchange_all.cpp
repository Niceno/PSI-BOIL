#include "vectorbool.h"

/******************************************************************************/
void VectorBool::exchange_all(const Comp comp, const int dir) {
	
  /* through all vector components */
  for_m(m)
    if(comp == Comp::undefined() || comp == m) 
      vec[m].exchange_all(dir);
}

/*-----------------------------------------------------------------------------+
 '$Id: vectorbool_exchange_all.cpp,v 1.1 2014/02/04 08:20:58 sato Exp $'/
+-----------------------------------------------------------------------------*/
