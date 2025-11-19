#include "nucleation.h"

/******************************************************************************/
bool Nucleation::height_bubble(const int ns ) {
/***************************************************************************//**
*  \brief calculate bottom of bubble height previously departed
*  true: bottom is higher than criterion zplant -> plant next bubble
*  false: bottom is lower than criterion zplant -> don't plant next bubble
*******************************************************************************/

  if(sites[ns].zplant()<0.0) {   // crude code
    return true;
  }

  /* bheight = true: bubble bottom is higher than zplant -> plant
   *           false: bubble bottom is lower than zplant -> don't plant */
  bool bheight;

  real zft=zftmin( Range<real>(sites[ns].x()-dxmin,sites[ns].x()+dxmin),
                   Range<real>(sites[ns].y()-dxmin,sites[ns].y()+dxmin),
                   Range<real>(-boil::exa,boil::exa));
  if (zft>sites[ns].zplant()) {
    bheight=true;
  } else {
    bheight=false;
  }

  return bheight;
}
