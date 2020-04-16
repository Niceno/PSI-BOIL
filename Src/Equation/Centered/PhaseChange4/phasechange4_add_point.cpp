#include "phasechange4.h"

/******************************************************************************/
void PhaseChange4::add_point(const int i0, const int j0, const int k0,
                             const int i1, const int j1, const int k1,
                             const Sign dir, const Comp & m,
                             const bool is_solid, bool & terminate,
                             std::vector<real> & stencil,
                             std::vector<real> & values) {
/***************************************************************************//**
*  \brief Add a point to a stencil.
*******************************************************************************/

  /* is there an interface? */
  if(!is_solid&&interface(dir,m,i0,j0,k0)) {

    terminate = true; 

    real tgamma;
    real dist_int = std::max(epsl*distance_center(dir,m,i0,j0,k0),
                            distance_int(dir,m,i0,j0,k0,tgamma));

    /* directional choice */
    if(dir>0) {
      stencil.push_back( dist_int);
    } else {
      stencil.push_back(-dist_int);
    }
    values.push_back(tgamma);

  /* are we at a solid-fluid boundary? */
  } else if(is_solid!=dom->ibody().off(i1,j1,k1)) {

    terminate = true; 
     
    /* directional choice */
    if(dir>0) {
      stencil.push_back( distance_face(dir,m,i0,j0,k0));
      values.push_back(bndtpr[m][i1][j1][k1]);
    } else {
      stencil.push_back(-distance_face(dir,m,i0,j0,k0));
      values.push_back(bndtpr[m][i0][j0][k0]);
    }

  /* neither interface, nor boundary */
  } else {
    
    /* directional choice */
    if(dir>0) {
      stencil.push_back( distance_center(dir,m,i0,j0,k0));
    } else {
      stencil.push_back(-distance_center(dir,m,i0,j0,k0));
    }
    values.push_back(tpr[i1][j1][k1]);

    /* are we at domain edge? */
    if(edge(dir,m,i0,j0,k0)) {
      terminate = true;
    }

  }

  return;
}

