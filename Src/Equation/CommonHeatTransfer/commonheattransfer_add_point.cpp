#include "commonheattransfer.h"

/******************************************************************************/
bool CommonHeatTransfer::add_point(const int i0, const int j0, const int k0,
                                   const int i1, const int j1, const int k1,
                                   const Sign dir, const Comp & m,
                                   const bool is_solid, bool & terminate,
                                   bool & interface_reached,
                                   std::vector<StencilPoint> & stencil,
                                   const Old old) const {
/***************************************************************************//**
*  \brief Add a point to a stencil. Output: if the previous point should be 
*         discarded since the interface is too close to it.
*******************************************************************************/

  /* is there an interface? */
  if(!is_solid&&interface(dir,m,i0,j0,k0,old)) {

    terminate = true; 
    interface_reached = true;

    real tgamma;
    real dist_int = distance_int(dir,m,i0,j0,k0,tgamma,old);

    StencilPoint sp(stencil.size(),tgamma,dist_int);
    /* directional choice */
    if(dir<0)
      sp.pos = -dist_int;
    stencil.push_back(sp);

    /* should the previous point be discarded? */
    if(dist_int<distance_face(dir,m,i0,j0,k0)) {
      return true;
    }

  /* are we at a solid-fluid boundary? */
  } else if(is_solid!=topo->domain()->ibody().off(i1,j1,k1)) {

    terminate = true; 

    StencilPoint sp(stencil.size());
     
    /* directional choice */
    if(dir>0) {
      sp.pos =  distance_face(dir,m,i0,j0,k0);
      if(is_solid)
        sp.val = bndtpr_sol[m][i1][j1][k1];
      else
        sp.val = bndtpr_flu[m][i1][j1][k1];
    } else {
      sp.pos = -distance_face(dir,m,i0,j0,k0);
      if(is_solid)
        sp.val = bndtpr_sol[m][i0][j0][k0];
      else
        sp.val = bndtpr_flu[m][i0][j0][k0];
    }
    stencil.push_back(sp);


  /* neither interface, nor boundary */
  } else {

    StencilPoint sp(stencil.size());
    
    /* directional choice */
    if(dir>0) {
      sp.pos =  distance_center(dir,m,i0,j0,k0);
    } else {
      sp.pos = -distance_center(dir,m,i0,j0,k0);
    }
    sp.val = tpr[i1][j1][k1];
    stencil.push_back(sp);

    /* are we at domain edge? */
    if(edge(dir,m,i0,j0,k0)) {
      terminate = true;
    }

  }

  return false;
}

