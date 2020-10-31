#include "commonheattransfer.h"

/******************************************************************************/
void CommonHeatTransfer::construct_stencil(
                                 std::vector<real> & stencil,
                                 std::vector<real> & values,
                                 const bool is_solid, const Comp & m,
                                 const int i, const int j, const int k,
                                 const AccuracyOrder & accuracy_order,
                                 const bool discard_points,
                                 const Old old) const {
/***************************************************************************//**
*  \brief Construct stencil for difference calculations.
*******************************************************************************/

  /* reset vectors */
  stencil.clear();
  values.clear();

  /* reference indices */
  int c_idx(-1), w_idx(-1), e_idx(-1), 
                 ww_idx(-1), ee_idx(-1),
                 www_idx(-1), eee_idx(-1);

  /* add zeroth point */
  stencil.push_back(0);
  values.push_back(tpr[i][j][k]);
  c_idx = stencil.size()-1;

  /* set-up generic directional indexing */
  int ii0(0),jj0(0),kk0(0);
  int ii1(0),jj1(0),kk1(0);

  int * idx0, * idx1;

  if       (m==Comp::i()) {
    idx0 = &ii0;
    idx1 = &ii1;
  } else if(m==Comp::j()) {
    idx0 = &jj0;
    idx1 = &jj1;
  } else {
    idx0 = &kk0;
    idx1 = &kk1;
  }

  /* set-up termination flags */
  bool wend(false), eend(false);
  bool interface_west(false), interface_east(false);
  int int_west_idx(-1), int_east_idx(-1);

  /* prepare possible discarding set */
  bool discard_center(false);
  std::set<int> discard_set = {};

  /*** add first and second point: always possible ***/
  
  /* west */
  *idx0 = 0;
  *idx1 = -1;

  bool discard_center_by_west = add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                                          Sign::neg(), m, is_solid, 
                                          wend, interface_west,
                                          stencil, values, old);
  w_idx = stencil.size()-1;

  /* interface position is recorded, if true */
  if(interface_west)
    int_west_idx = w_idx;

  /* east */
  *idx0 = 0;
  *idx1 = +1;

  bool discard_center_by_east = add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                                          Sign::pos(), m, is_solid, 
                                          eend, interface_east,
                                          stencil, values, old);
  e_idx = stencil.size()-1;

  /* interface position is recorded, if true */
  if(interface_east)
    int_east_idx = e_idx;

  /* special treatment for first order */
  if(accuracy_order.eval()==1) {
    int dscrd_idx(-1);
    /* 1. upwind AND interface on only one side */
    if(accuracy_order.upwind()&&(interface_east!=interface_west)) {
      /* if interface in west, discard east and vice versa */
      dscrd_idx = interface_west*e_idx+interface_east*w_idx;
    /* 2. central  */
    } else {
      dscrd_idx = c_idx;
    }

    stencil.erase(stencil.begin() + dscrd_idx);
    values.erase(values.begin() + dscrd_idx);
    return;

  } else {
    /* either marks center for discarding */
    discard_center = discard_center_by_west | discard_center_by_east;
  }

  /* check for discarding */
  if(discard_center)
    discard_set.insert(c_idx);

  /*** if termination flags did not fire, add further points ***/

  /* west-west */
  if(!wend) {
    *idx0 = -1;
    *idx1 = -2;
  
    /* w-w point possibly marks w point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::neg(), m, is_solid, wend, interface_west,
                 stencil, values, old))
      discard_set.insert(w_idx);

    ww_idx = stencil.size()-1;

#if 0
    /* interface position is recorded, if true */
    if(interface_west)
      int_west_idx = ww_idx;
#endif

    /* correct distance */
    stencil[ww_idx] += stencil[w_idx];
  }

  /* east-east */
  if(!eend) {
    *idx0 = +1;
    *idx1 = +2;
  
    /* e-e point possibly marks e point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::pos(), m, is_solid, eend, interface_east,
                 stencil, values, old))
      discard_set.insert(e_idx);

    ee_idx = stencil.size()-1;

#if 0
    /* interface position is recorded, if true */
    if(interface_east)
      int_east_idx = ee_idx;
#endif

    /* correct distance */
    stencil[ee_idx] += stencil[e_idx];
  }

  /*** if termination flags did not fire, add further points ***/

  /* west-west-west */
  if(!wend) {
    *idx0 = -2;
    *idx1 = -3;

    /* w-w-w point possibly marks w-w point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::neg(), m, is_solid, wend, interface_west,
                 stencil, values, old))
      discard_set.insert(ww_idx);

    www_idx = stencil.size()-1;

#if 0
    /* interface position is recorded, if true */
    if(interface_west)
      int_west_idx = www_idx;
#endif

    /* correct distance */
    stencil[www_idx] += stencil[ww_idx];
  }

  /* east-east-east */
  if(!eend) {
    *idx0 = +2;
    *idx1 = +3;

    /* e-e-e point possibly marks e-e point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::pos(), m, is_solid, eend, interface_east,
                 stencil, values, old))
      discard_set.insert(ee_idx);

    eee_idx = stencil.size()-1;

#if 0
    /* interface position is recorded, if true */
    if(interface_east)
      int_east_idx = eee_idx;
#endif

    /* correct distance */
    stencil[eee_idx] += stencil[ee_idx];
  }

  /*** do we want to reposition the stencil to an interface? ***/
  /* only used if upwind_to_interface is true and interface is directly
     adjacent from one side -> xor gate used */
  if(accuracy_order.upwind()&&(interface_west^interface_east)) {
    real new_origin = interface_west ?
                      stencil[int_west_idx] :
                      stencil[int_east_idx];
    
    for(int idx(0); idx != stencil.size(); ++idx) {
      stencil[idx] -= new_origin;
    }
  }

  /*** the stencil is now 3-7 points, ordered by importance ***/
  
  /*** do we want to discard points too close to an interface? ***/
  if(discard_points) {
    /* A set is in ascending order as per STL and we can remove elements by
       using reverse iteration */
    for(auto rit = discard_set.rbegin(); rit != discard_set.rend(); rit++) {
      stencil.erase( stencil.begin() + (*rit) );
      values.erase( values.begin() + (*rit) );
    }
  }

  return;
}
