#include "additive.h"

/******************************************************************************/
void AC::v_cycle(const int l) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_solver(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    v_cycle(l+1);
 
    /* post-smooth */
    call_solver(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));

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
void AC::f_cycle(const int l, const bool upstream) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_solver(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    f_cycle(l+1,upstream);

    /* re-smooth */
    call_solver(l,MaxIter(10 * (l+1)),ResRat(0.01),ResTol(boil::pico));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /*
       prolongate residual between the two recursive calls
       (not in original F-cycle)
    */
    if(upstream && l > 0)
      interpolation(*L[l], *L[l-1]);
    
    /* recursive call to a coarser level on the way back */
    v_cycle(l+1);

    /* post-smooth */
    call_solver(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));

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
void AC::w_cycle(const int l, const bool upstream) {

  if(l!=nlevels-1) {
    /* pre-smooth */
    call_solver(l,MaxIter(5),ResRat(0.1),ResTol(boil::nano));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    w_cycle(l+1,upstream);

    /* re-smooth */
    call_solver(l,MaxIter(10 * (l+1)),ResRat(0.01),ResTol(boil::pico));

    /* restrict residual to coarser grid */
    residual(*L[l]);
    restriction(*L[l], *L[l+1]);

    /*
       prolongate residual between the two recursive calls
       (not in original W-cycle)
    */
    if(upstream && l > 0)
      interpolation(*L[l], *L[l-1]);
    
    /* recursive call to a coarser level on the way back */
    w_cycle(l+1,upstream);

    /* post-smooth */
    call_solver(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));

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
void AC::full_cycle(const int l, const Cycle & cyc) {

  if(l!=nlevels-1) {
    /* restrict fnew to coarser grid */
    restriction(*L[l], *L[l+1]);

    /* recursive call to a coarser level */
    full_cycle(l+1,cyc);

    /* recursive call to a coarser level on the way back */
    if       (cyc==Cycle::V()) {
      v_cycle(l);
    } else if(cyc==Cycle::F1()) {
      f_cycle(l,false);
    } else if(cyc==Cycle::F2()) {
      f_cycle(l,true);
    } else if(cyc==Cycle::W1()) {
      w_cycle(l,false);
    } else if(cyc==Cycle::W2()) {
      w_cycle(l,true);
      /* just solve (shouldn't happen) */
    } else {
      call_solver(l,MaxIter(20 * (l+1)),ResRat(0.01),ResTol(boil::pico));
    }

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

