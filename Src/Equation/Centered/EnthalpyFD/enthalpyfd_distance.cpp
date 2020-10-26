#include "enthalpyfd.h"

/******************************************************************************/
real EnthalpyFD::distance_center(const Sign sig, const Comp & m,
                                 const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell center  
*******************************************************************************/
  if       (m==Comp::i()) {
    if(sig<0) {
      return phi.dxw(i);
    } else {
      return phi.dxe(i);
    }
  } else if(m==Comp::j()) {
    if(sig<0) {
      return phi.dys(j);
    } else {
      return phi.dyn(j);
    }
  } else {
    if(sig<0) {
      return phi.dzb(k);
    } else {
      return phi.dzt(k);
    }
  }

  return 0.0;
}

/******************************************************************************/
real EnthalpyFD::distance_face(const Sign sig, const Comp & m,
                               const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell face  
*******************************************************************************/
  if       (m==Comp::i()) {
    return 0.5*phi.dxc(i);
  } else if(m==Comp::j()) {
    return 0.5*phi.dyc(j);
  } else {
    return 0.5*phi.dzc(k);
  }

  return 0.0;
}

/***************************************************************************//*** 
*  \brief calculate distance to interface, as well as interface tpr  
*******************************************************************************/
/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
*/

/*** new ***/

/* generic */
real EnthalpyFD::distance_int(const Sign dir, const Comp & m,
                              const int i, const int j, const int k,
                              real & tint) {
  if        (m==Comp::i()) {
    return distance_int_x(dir,i,j,k,tint);
  } else if (m==Comp::j()) {
    return distance_int_y(dir,i,j,k,tint);
  } else {
    return distance_int_z(dir,i,j,k,tint);
  }

  return 0.0;
}


/* x-direction */
real EnthalpyFD::distance_int_x(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_x(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint(i,j,k);
    dist += ghost_distance(Comp::i(),cell_marker,i,j,k);
  } else {
    tint = Tint(i+int(dir),j,k);
    dist += ghost_distance(Comp::i(),cell_marker,i+int(dir),j,k);
  }

  return dist;
}

/* y-direction */
real EnthalpyFD::distance_int_y(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_y(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint(i,j,k);
    dist += ghost_distance(Comp::j(),cell_marker,i,j,k);
  } else {
    tint = Tint(i,j+int(dir),k);
    dist += ghost_distance(Comp::j(),cell_marker,i,j+int(dir),k);
  }

  return dist;
}

/* z-direction */
real EnthalpyFD::distance_int_z(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_z(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint(i,j,k);
    dist += ghost_distance(Comp::k(),cell_marker,i,j,k);
  } else {
    tint = Tint(i,j,k+int(dir));
    dist += ghost_distance(Comp::k(),cell_marker,i,j,k+int(dir));
  }

  return dist;
}

/*** old ***/

/* generic */
real EnthalpyFD::distance_int_old(const Sign dir, const Comp & m,
                              const int i, const int j, const int k,
                              real & tint) {
  if        (m==Comp::i()) {
    return distance_int_x_old(dir,i,j,k,tint);
  } else if (m==Comp::j()) {
    return distance_int_y_old(dir,i,j,k,tint);
  } else {
    return distance_int_z_old(dir,i,j,k,tint);
  }

  return 0.0;
}


/* x-direction */
real EnthalpyFD::distance_int_x_old(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_x_old(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint_old(i,j,k);
    dist += ghost_distance(Comp::i(),cell_marker,i,j,k);
  } else {
    tint = Tint_old(i+int(dir),j,k);
    dist += ghost_distance(Comp::i(),cell_marker,i+int(dir),j,k);
  }

  return dist;
}

/* y-direction */
real EnthalpyFD::distance_int_y_old(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_y_old(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint_old(i,j,k);
    dist += ghost_distance(Comp::j(),cell_marker,i,j,k);
  } else {
    tint = Tint_old(i,j+int(dir),k);
    dist += ghost_distance(Comp::j(),cell_marker,i,j+int(dir),k);
  }

  return dist;
}

/* z-direction */
real EnthalpyFD::distance_int_z_old(const Sign dir, 
                                const int i, const int j, const int k,
                                real & tint) {
  Sign cell_marker;
  real dist = topo->distance_int_z_old(dir,i,j,k,cell_marker);
  if(cell_marker < 0) {
    tint = Tint_old(i,j,k);
    dist += ghost_distance(Comp::k(),cell_marker,i,j,k);
  } else {
    tint = Tint_old(i,j,k+int(dir));
    dist += ghost_distance(Comp::k(),cell_marker,i,j,k+int(dir));
  }

  return dist;
}
