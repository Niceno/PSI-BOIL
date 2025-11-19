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
