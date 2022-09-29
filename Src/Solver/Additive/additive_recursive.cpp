#include "additive.h"

/******************************************************************************/
void AC::v_cycle(const int l,const std::array<MaxIter,3> & mv,
                 const std::array<MaxIter,3> & ms, const int gi) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_smoother(l,MinIter(1),mv[0],resrat_val,restol_val,ms[0],gi);

    /* restrict residual to coarser grid */
    //residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    v_cycle(l+1,mv,ms,gi);
 
    /* post-smooth */
    call_solver(l,MinIter(1),mv[2],resrat_val,restol_val,ms[2],gi);

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::f_cycle(const int l,const std::array<MaxIter,3> & mv,
                 const std::array<MaxIter,3> & ms, const int gi) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_smoother(l,MinIter(1),mv[0],resrat_val,restol_val,ms[0],gi);

    /* restrict residual to coarser grid */
    //residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    f_cycle(l+1,mv,ms,gi);

    /* re-smooth */
    call_smoother(l,MinIter(1),mv[1],resrat_val,restol_val,ms[1],gi);

    /* restrict residual to coarser grid */
    //residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level on the way back */
    v_cycle(l+1,mv,ms,gi);

    /* post-smooth */
    call_solver(l,MinIter(1),mv[2],resrat_val,restol_val,ms[2],gi);

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::w_cycle(const int l,const std::array<MaxIter,3> & mv,
                 const std::array<MaxIter,3> & ms, const int gi) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_smoother(l,MinIter(1),mv[0],resrat_val,restol_val,ms[0],gi);

    /* restrict residual to coarser grid */
    //residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    w_cycle(l+1,mv,ms,gi);

    /* re-smooth */
    call_smoother(l,MinIter(1),mv[1],resrat_val,restol_val,ms[1],gi);

    /* restrict residual to coarser grid */
    //residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level on the way back */
    w_cycle(l+1,mv,ms,gi);

    /* post-smooth */
    call_solver(l,MinIter(1),mv[2],resrat_val,restol_val,ms[2],gi);

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::full_cycle(const int l, const Cycle & cyc,
                    const std::array<MaxIter,3> & mv,
                    const std::array<MaxIter,3> & ms,
                    const int gi) {

  if(l!=nlevels-1) {
    /* restrict fnew to coarser grid */
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    full_cycle(l+1,cyc,mv,ms,gi);

    /* recursive call to a coarser level on the way back */
    if       (cyc==Cycle::V()) {
      v_cycle(l,mv,ms,gi);
    } else if(cyc==Cycle::F()) {
      f_cycle(l,mv,ms,gi);
    } else if(cyc==Cycle::W()) {
      w_cycle(l,mv,ms,gi);
    } else if(cyc==Cycle::flex()) {
      flex_cycle(l,mv,ms,gi);
      /* just solve (shouldn't happen) */
    } else {
      call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);
    }

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::flex_cycle(const int l,const std::array<MaxIter,3> & mv,
                    const std::array<MaxIter,3> & ms, 
                    const int gi, const Sign sig) {

  int idx = 0;
  if(sig == Sign::pos())
    idx = 2;

  if(l!=nlevels-1) {
    /* smooth until convergence or stale */
    if(call_smoother(l,MinIter(1),mv[idx],resrat_val,restol_val,ms[idx],gi)) {
      if(l==0) {
        return;
      } else {
        /* prolongate to a finer level */
        interpolation(*L[l], *L[l-1]);

        /* recursive call to a finer level */
        flex_cycle(l-1,mv,ms,gi,Sign::pos());
      }
    } else {
      /* restrict residual to coarser grid */
      restriction(*L[l], *L[l+1]);

      /* recursive call to a coarser level */
      flex_cycle(l+1,mv,ms,gi,Sign::neg());
    }
  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MinIter(1),MaxIter(100),ResRat(-1),restol_val,MaxIter(-1),gi);

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
    
    /* recursive call to a finer level */
    flex_cycle(l-1,mv,ms,gi);
  }

  return;
}
