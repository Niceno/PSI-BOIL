#include "profile.h"

/******************************************************************************/
void Profile::scale_coords(const real & s) {

  for(int i=0; i<coord.size(); i++) 
    coord[i] *= s;   

}

/******************************************************************************/
void Profile::scale_values(const real & s) {

  for(int i=0; i<value.size(); i++) 
    value[i] *= s;   

}
