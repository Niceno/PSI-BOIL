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

    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

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

    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

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

    if(l > 0)
      interpolation(*L[l], *L[l-1]);

  } else {
    /* solve 'precisely' at the coarsest level */
    call_solver(l,MaxIter(40 * (l+1)),ResRat(0.001),ResTol(boil::femto));

    interpolation(*L[l], *L[l-1]);
  }

  return;
}
