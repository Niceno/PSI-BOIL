#include "property.h"

/*============================================================================*/
void Property::varies(const Domain & d) {

  con = false;

  dom = & d;
	
  allocate(dom->ni(), dom->nj(), dom->nk());

  bndcnd = new BndCnd( *dom ); 

  for_aijk(i,j,k)
    val[i][j][k] = cval;

}

/*-----------------------------------------------------------------------------+
 '$Id: property_varies.cpp,v 1.2 2011/05/29 09:45:37 niceno Exp $'/
+-----------------------------------------------------------------------------*/
