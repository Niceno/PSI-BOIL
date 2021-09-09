#include "commonheattransfer.h"

/******************************************************************************/
real CommonHeatTransfer::distance_center(const Sign sig, const Comp & m,
                                         const int i, const int j, const int k)
                                         const {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell center  
*******************************************************************************/
  if       (m==Comp::i()) {
    if(sig<0) {
      return topo->clrold.dxw(i);
    } else {
      return topo->clrold.dxe(i);
    }
  } else if(m==Comp::j()) {
    if(sig<0) {
      return topo->clrold.dys(j);
    } else {
      return topo->clrold.dyn(j);
    }
  } else {
    if(sig<0) {
      return topo->clrold.dzb(k);
    } else {
      return topo->clrold.dzt(k);
    }
  }

  return 0.0;
}

/******************************************************************************/
real CommonHeatTransfer::distance_face(const Sign sig, const Comp & m,
                                       const int i, const int j, const int k)
                                       const {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell face  
*******************************************************************************/
  if       (m==Comp::i()) {
    return 0.5*topo->clrold.dxc(i);
  } else if(m==Comp::j()) {
    return 0.5*topo->clrold.dyc(j);
  } else {
    return 0.5*topo->clrold.dzc(k);
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
real CommonHeatTransfer::distance_int(const Sign dir, const Comp & m,
                                      const int i, const int j, const int k,
                                      real & tint, Sign & cell_marker,
                                      const ResistEval re,
                                      const Old old) const {
  if        (m==Comp::i()) {
    return distance_int_x(dir,i,j,k,tint,cell_marker,re,old);
  } else if (m==Comp::j()) {
    return distance_int_y(dir,i,j,k,tint,cell_marker,re,old);
  } else {
    return distance_int_z(dir,i,j,k,tint,cell_marker,re,old);
  }

  return 0.0;
}


/* x-direction */
real CommonHeatTransfer::distance_int_x(const Sign dir, 
                                        const int i, const int j, const int k,
                                        real & tint, Sign & cell_marker, 
                                        const ResistEval re,
                                        const Old old) const {
  real dist = topo->distance_int_x(dir,i,j,k,cell_marker,old);
  if(cell_marker < 0) {
    tint = Tint(i,j,k,old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::i(),i,j,k,  i,j,k,   i-int(dir),j,k,
              dist,cell_marker,tint,old); 
    }
  } else {
    tint = Tint(i+int(dir),j,k,old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::i(),i,j,k,  i+int(dir),j,k,   i-int(dir),j,k,
              dist,cell_marker,tint,old); 
    }
  }

  return dist;
}

/* y-direction */
real CommonHeatTransfer::distance_int_y(const Sign dir, 
                                        const int i, const int j, const int k,
                                        real & tint, Sign & cell_marker, 
                                        const ResistEval re,
                                        const Old old) const {
  real dist = topo->distance_int_y(dir,i,j,k,cell_marker,old);
  if(cell_marker < 0) {
    tint = Tint(i,j,k,old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::j(),i,j,k,   i,j,k,   i,j-int(dir),k,
              dist,cell_marker,tint,old); 
    }
  } else {
    tint = Tint(i,j+int(dir),k,old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::j(),i,j,k,   i,j+int(dir),k,   i,j-int(dir),k,
              dist,cell_marker,tint,old); 
    }
  }

  return dist;
}

/* z-direction */
real CommonHeatTransfer::distance_int_z(const Sign dir, 
                                        const int i, const int j, const int k,
                                        real & tint, Sign & cell_marker, 
                                        const ResistEval re,
                                        const Old old) const {
  real dist = topo->distance_int_z(dir,i,j,k,cell_marker,old);
  if(cell_marker < 0) {
    tint = Tint(i,j,k,old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::k(),i,j,k,   i,j,k,   i,j,k-int(dir),
              dist,cell_marker,tint,old); 
    }
  } else {
    tint = Tint(i,j,k+int(dir),old);

    /* interfacial heat transfer resistance */
    if(use_int_resist&&re==ResistEval::yes) {
      resTint(-dir,Comp::k(),i,j,k,   i,j,k+int(dir),   i,j,k-int(dir),
              dist,cell_marker,tint,old); 
    }
  }

  return dist;
}
