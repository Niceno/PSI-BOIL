#include "simulation_time.h"

/******************************************************************************/
Times::Times() {
/*-----------------------------------------------------------------------------+ 
|  nothing is specified -> steady                                              |
+-----------------------------------------------------------------------------*/
  f_dt   = 0;
  t_dt   = 1; 
  d_t    = 1.0e+30; 
  t_time = d_t * t_dt;
  d_ti   = 1.0/d_t; 
  time   = 0.0;
  nhist  = 0;
  istep_inc = 0;
  coef_inc=1.2;
  coef_dec=0.4;
  print_cinfo = true;
  alloc1d(& cfl_hist, CFL_HIST_SIZE );

  print_total();
} 

/******************************************************************************/
Times::Times(const int n, const real d) {
/*-----------------------------------------------------------------------------+ 
|  number of time steps and time step are specified                            |
+-----------------------------------------------------------------------------*/
  f_dt   = 0;
  t_dt   = n; 
  d_t    = d; 
  t_time = d_t * t_dt;
  d_ti   = 1.0/d_t;
  time   = 0.0;
  nhist  = 0;
  istep_inc = 0;
  coef_inc=1.2;
  coef_dec=0.4;
  print_cinfo = true;
  alloc1d(& cfl_hist, CFL_HIST_SIZE );

  print_total();
} 

/******************************************************************************/
Times::Times(const real t, const real d) {
/*-----------------------------------------------------------------------------+ 
|  total time and time step are specified                                      |
+-----------------------------------------------------------------------------*/
  t_time = t;
  f_dt   = 0; 
  d_t    = d; 
  t_dt   = (int)(t/d); 
  d_ti   = 1.0/d_t;
  time   = 0.0;
  nhist  = 0;
  istep_inc = 0;
  coef_inc=1.2;
  coef_dec=0.4;
  print_cinfo = true;
  alloc1d(& cfl_hist, CFL_HIST_SIZE );

  print_total();
} 

/******************************************************************************/
Times::Times(const real t, const int n) {
/*-----------------------------------------------------------------------------+ 
|  total time and number of time steps are specified                           |
+-----------------------------------------------------------------------------*/
  t_time = t;
  f_dt   = 0; 
  d_t    = t/(real)n; 
  t_dt   = n; 
  d_ti   = 1.0/d_t;
  time   = 0.0;
  nhist  = 0;
  istep_inc = 0;
  coef_inc=1.2;
  coef_dec=0.4;
  print_cinfo = true;
  alloc1d(& cfl_hist, CFL_HIST_SIZE );

  print_total();
} 

/*-----------------------------------------------------------------------------+
 '$Id: simulation_time.cpp,v 1.5 2016/05/05 08:04:35 sato Exp $'/
+-----------------------------------------------------------------------------*/
