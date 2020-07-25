#include "cavitypressure.h"

/***************************************************************************//*** 
*  \brief calculate distance to interface, as well as interface pressure  
*******************************************************************************/
/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
*/

/* generic */
real CavityPressure::distance_int(const Sign dir, const Comp & m,
                                  const int i, const int j, const int k,
                                  real & pint) {
  if        (m==Comp::i()) {
    return distance_int_x(dir,i,j,k,pint);
  } else if (m==Comp::j()) {
    return distance_int_y(dir,i,j,k,pint);
  } else {
    return distance_int_z(dir,i,j,k,pint);
  }

  return 0.0;
}


/* x-direction */
real CavityPressure::distance_int_x(const Sign dir, 
                                    const int i, const int j, const int k,
                                    real & pint) {
  Sign cell_marker;
  real dist = topo->distance_int_x(dir,i,j,k,cell_marker);
#if 1
  if(cell_marker < 0) {
    pint = Pint_wrapper(i,j,k);
  } else {
    pint = Pint_wrapper(i+int(dir),j,k);
  }
#else
  pint = 0.5*(Pint_wrapper(i,j,k)+Pint_wrapper(i+int(dir),j,k));
#endif

  return dist;
}

/* y-direction */
real CavityPressure::distance_int_y(const Sign dir, 
                                    const int i, const int j, const int k,
                                    real & pint) {
  Sign cell_marker;
  real dist = topo->distance_int_y(dir,i,j,k,cell_marker);
#if 1
  if(cell_marker < 0) {
    pint = Pint_wrapper(i,j,k);
  } else {
    pint = Pint_wrapper(i,j+int(dir),k);
  }
#else
  pint = 0.5*(Pint_wrapper(i,j,k)+Pint_wrapper(i,j+int(dir),k));
#endif

  return dist;
}

/* z-direction */
real CavityPressure::distance_int_z(const Sign dir, 
                                    const int i, const int j, const int k,
                                    real & pint) {
  Sign cell_marker;
  real dist = topo->distance_int_z(dir,i,j,k,cell_marker);
#if 1
  if(cell_marker < 0) {
    pint = Pint_wrapper(i,j,k);
  } else {
    pint = Pint_wrapper(i,j,k+int(dir));
  }
#else
  pint = 0.5*(Pint_wrapper(i,j,k)+Pint_wrapper(i,j,k+int(dir)));
#endif

  return dist;
}
