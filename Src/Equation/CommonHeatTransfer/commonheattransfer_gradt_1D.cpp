#include "commonheattransfer.h"

/******************************************************************************/
real CommonHeatTransfer::gradt1D(const bool is_solid, const Comp & m,
                                 const int i, const int j, const int k,
                                 const AccuracyOrder & accuracy_order,
                                 const bool discard_points) const {
/***************************************************************************//**
*  \brief Calculate grad(tpr) in a given cell.
*******************************************************************************/

  /* set-up arrays */
  std::vector<real> stencil;
  std::vector<real> values;

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

  /* prepare possible discarding set */
  bool discard_center(false);
  std::set<int> discard_set = {};

  /* for first order schemes, center is always discarded */
  if(accuracy_order.eval()<2) {
    discard_center = true;
  }

  /*** add first and second point: always possible ***/
  
  /* west */
  *idx0 = 0;
  *idx1 = -1;

  discard_center |= add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                              Sign::neg(), m, is_solid, wend, stencil, values);
  w_idx = stencil.size()-1;

  /* east */
  *idx0 = 0;
  *idx1 = +1;

  discard_center |= add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                              Sign::pos(), m, is_solid, eend, stencil, values);
  e_idx = stencil.size()-1;

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
                 Sign::neg(), m, is_solid, wend, stencil, values))
      discard_set.insert(w_idx);

    ww_idx = stencil.size()-1;

    /* correct distance */
    stencil[ww_idx] += stencil[w_idx];
  }

  /* east-east */
  if(!eend) {
    *idx0 = +1;
    *idx1 = +2;
  
    /* e-e point possibly marks e point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::pos(), m, is_solid, eend, stencil, values))
      discard_set.insert(e_idx);

    ee_idx = stencil.size()-1;

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
                 Sign::neg(), m, is_solid, wend, stencil, values))
      discard_set.insert(ww_idx);

    www_idx = stencil.size()-1;

    /* correct distance */
    stencil[www_idx] += stencil[ww_idx];
  }

  /* east-east-east */
  if(!eend) {
    *idx0 = +2;
    *idx1 = +3;

    /* e-e-e point possibly marks e-e point for discard */
    if(add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                 Sign::pos(), m, is_solid, eend, stencil, values))
      discard_set.insert(ee_idx);

    eee_idx = stencil.size()-1;

    /* correct distance */
    stencil[eee_idx] += stencil[ee_idx];
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

  /*** the stencil is now 2-7 points, ordered by importance ***/
  int diff_req = stencil.size()-1;
  AccuracyOrder ao(std::min(diff_req,accuracy_order.eval()));

  /* differences are implemented up to fourth-order */
#if 0
  //if(m==Comp::i()&& (iflag[i][j][k]==1||iflag[i][j][k]==2) ) {
  if( (i==5&&k==28) || (i==28&&k==5) || (i==10&&k==27) || (i==27&&k==10) ) {
    boil::oout<<"gradt1D: "<<i<<" "<<k<<" "<<m<<" "<<iflag[i][j][k]<<" |";
    for(auto s : stencil) 
      boil::oout<<" "<<s;
    boil::oout<<" |";
    for(auto v : values) 
      boil::oout<<" "<<v;
    boil::oout<<" | "
              <<topo->nth_order_difference(stencil,values,std::min(accuracy_order,diff_req))
              <<boil::endl;
  }
#endif

  return topo->nth_order_difference(stencil,values,ao);

}
