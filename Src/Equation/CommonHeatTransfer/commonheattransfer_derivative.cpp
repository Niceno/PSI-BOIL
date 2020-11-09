#include "commonheattransfer.h"

/******************************************************************************/
real CommonHeatTransfer::first_derivative(const bool is_solid, const Comp & m,
                                          const int i, const int j, const int k,
                                          const AccuracyOrder & accuracy_order,
                                          const bool discard_points,
                                          const Old old) const {
/***************************************************************************//**
*  \brief Calculate grad(tpr) in a given cell.
*******************************************************************************/

  /* set-up arrays */
  std::vector<StencilPoint> stencil;

  /* construct the stencil and fill it with values */
  construct_stencil(stencil,is_solid,m,i,j,k,
                    accuracy_order,discard_points,old);

  /*** the stencil is now 2-7 points, ordered by importance ***/
  int diff_req = stencil.size()-1;
  AccuracyOrder ao(std::min(diff_req,accuracy_order.eval()));

  /* differences are implemented up to fourth-order */

  return topo->nth_order_first(stencil,ao);

}

/******************************************************************************/
real CommonHeatTransfer::second_derivative(const bool is_solid, const Comp & m,
                                          const int i, const int j, const int k,
                                          const AccuracyOrder & accuracy_order,
                                          const bool discard_points,
                                          const Old old) const {
/***************************************************************************//**
*  \brief Calculate grad(tpr) in a given cell.
*******************************************************************************/

  /* set-up arrays */
  std::vector<StencilPoint> stencil;

  /* construct the stencil and fill it with values */
  construct_stencil(stencil,is_solid,m,i,j,k,
                    accuracy_order,discard_points,old);

  /*** the stencil is now 2-7 points, ordered by importance ***/
  int diff_req = stencil.size()-1;
  AccuracyOrder ao(std::min(diff_req,accuracy_order.eval()));

  /* differences are implemented up to fourth-order */

  return topo->nth_order_second(stencil,ao);

}
