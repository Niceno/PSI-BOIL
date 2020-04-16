#include "phasechange4.h"

/******************************************************************************/
real PhaseChange4::gradt1D(const bool is_solid, const Comp & m,
                           const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief Calculate grad(tpr) in a given cell.
*******************************************************************************/

  /* set-up arrays */
  std::vector<real> stencil;
  std::vector<real> values;

  /* reference indices */
  int c_idx(-1), w_idx(-1), e_idx(-1), ww_idx(-1), ee_idx(-1);

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

  /*** add first and second point: always possible ***/
  
  /* west */
  *idx0 = 0;
  *idx1 = -1;

  add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
            Sign::neg(), m, is_solid, wend, stencil, values);
  w_idx = stencil.size()-1;

  /* east */
  *idx0 = 0;
  *idx1 = +1;

  add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
            Sign::pos(), m, is_solid, eend, stencil, values);
  e_idx = stencil.size()-1;

  /*** terminate if second-order accuracy is desired ***/
  if(use_second_order_accuracy)
    return second_order_difference(stencil,values);
   
  /*** if termination flags did not fire, add further points ***/

  /* west-west */
  if(!wend) {
    *idx0 = -1;
    *idx1 = -2;
  
    add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
              Sign::neg(), m, is_solid, wend, stencil, values);
    ww_idx = stencil.size()-1;

    /* correct distance */
    stencil[ww_idx] += stencil[w_idx];
  }

  /* east-east */
  if(!eend) {
    *idx0 = +1;
    *idx1 = +2;
  
    add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
              Sign::pos(), m, is_solid, eend, stencil, values);
    ee_idx = stencil.size()-1;

    /* correct distance */
    stencil[ee_idx] += stencil[e_idx];
  }

  /*** check the stencil size ***/

  /* full stencil = 5 points */
  if(stencil.size()==5) {
    return fourth_order_difference(stencil,values);

  /* minimal stencil = 3 points */
  } else if(stencil.size()==3) {
    return second_order_difference(stencil,values);

  /* intermediate value = 4 points */
  } else if(stencil.size()==4) {

    /* can we extend to 5? */
    if(wend&&eend) {
      /* no */
      return third_order_difference(stencil,values);

    /* yes */
    } else {
  
      /* catastrophic error */
      assert(wend!=eend);

      /* we can go further east */
      if(wend) {

        /* east-east-east */
        *idx0 = +2;
        *idx1 = +3;

        add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                  Sign::pos(), m, is_solid, eend, stencil, values);

        /* correct distance */
        stencil.back() += stencil[ee_idx];

      /* we can go further west */
      } else {

        /* west-west-west */
        *idx0 = -2;
        *idx1 = -3;

        add_point(i+ii0,j+jj0,k+kk0, i+ii1,j+jj1,k+kk1,
                  Sign::neg(), m, is_solid, wend, stencil, values);

        /* correct distance */
        stencil.back() += stencil[ww_idx];

      }

      return fourth_order_difference(stencil,values);

    }

  /* catastrophic error */
  } else {
    boil::aout<<"PC4:gradt1D Inconsistent stencil size! Exiting."
              <<boil::endl;
    exit(0);
  }

  return 0.0;
}

