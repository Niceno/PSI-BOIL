#include "additive.h"

/******************************************************************************/
void AC::v_cycle(const int l,const std::array<MaxIter,3> & mv) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    //call_smoother(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));
    call_smoother(l,mv[0],ResRat(boil::atto),ResTol(boil::atto));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    v_cycle(l+1,mv);
 
    /* post-smooth */
    //call_smoother(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    call_smoother(l,mv[2],ResRat(boil::atto),ResTol(boil::atto));

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(100),ResRat(1e-6),ResTol(boil::atto));

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::f_cycle(const int l,const std::array<MaxIter,3> & mv) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    //call_smoother(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));
    call_smoother(l,mv[0],ResRat(boil::atto),ResTol(boil::atto));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    f_cycle(l+1,mv);

    /* re-smooth */
    //call_smoother(l,MaxIter(10 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    call_smoother(l,mv[1],ResRat(boil::atto),ResTol(boil::atto));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level on the way back */
    v_cycle(l+1,mv);

    /* post-smooth */
    //call_smoother(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    call_smoother(l,mv[2],ResRat(boil::atto),ResTol(boil::atto));

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::w_cycle(const int l,const std::array<MaxIter,3> & mv) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    //call_smoother(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));
    call_smoother(l,mv[0],ResRat(boil::atto),ResTol(boil::atto));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    w_cycle(l+1,mv);

    /* re-smooth */
    //call_smoother(l,MaxIter(10 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    call_smoother(l,mv[1],ResRat(boil::atto),ResTol(boil::atto));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level on the way back */
    w_cycle(l+1,mv);

    /* post-smooth */
    //call_smoother(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    call_smoother(l,mv[2],ResRat(boil::atto),ResTol(boil::atto));

    /* prolongate to a finer level */
    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

/******************************************************************************/
void AC::full_cycle(const int l, const Cycle & cyc,
                    const std::array<MaxIter,3> & mv) {

  if(l!=nlevels-1) {
    /* restrict fnew to coarser grid */
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    full_cycle(l+1,cyc,mv);

    /* recursive call to a coarser level on the way back */
    if       (cyc==Cycle::V()) {
      v_cycle(l,mv);
    } else if(cyc==Cycle::F()) {
      f_cycle(l,mv);
    } else if(cyc==Cycle::W()) {
      w_cycle(l,mv);
      /* just solve (shouldn't happen) */
    } else {
      call_solver(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    }

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

    /* prolongate to a finer level */
    interpolation(*L[l], *L[l-1]);
  }

  return;
}

