/******************************************************************************/
/* ------------ boundary or initial conditions */

  /* temperatures */
  const real tsat0 = 0.0;
  const real twall = tsat0 + deltat_wall;
  const real tout  = tsat0 + deltat_out;
  const real tnucl = tsat0 + deltat_nucl;
  const real tsat0_K = IF97::Tsat97(prs);

  /* heater power */
  real qsrc;
  if(LZheat>0.) {
    qsrc = qflux/fabs(LZheat);  /* [W/m3] */
    boil::oout<<"qsrc (explicit)= "<<qsrc<<" "<<qflux<<"\n";
  } else {
    qsrc = 0.0; /* W/m2 */
    boil::oout<<"qsrc (dirac)= "<<qsrc<<" "<<qflux<<"\n";
  }

  /* other */
  const real gravity = boil::g;

/******************************************************************************/
/* ------------ numerical simulation settings */

  /* total number of steps */
  const int ndt = 10e6; /* inconsequential */
  
  /* plot every */
  real t_per_plot = 0.01*1e-3;
  if(case_flag==0)
    t_per_plot = 0.1;

  /* steps per backup */
  const int n_per_backup = 5000;

  /* dt settings */
  const real surftens_dt_coef = 10.;
  const real initdtcoef = 1./10.;

  /* cfl with and without interfaces */
  real cfl_limit(0.1);
  if(case_flag==0)
    cfl_limit = 0.2;

  /* only liquid beyond this (fractional height) */
  const real zmax_mult = 0.9;

/* ------------ optional simulation settings */

  /* multigrid */
  const bool multigrid_stop_if_diverging = true;
  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 20;

  const int multigrid_niter = 30;
  MaxIter multigrid_mm = MaxIter(multigrid_niter);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm,multigrid_mm,multigrid_mm};

  ResRat multigrid_rr = ResRat(5e-5);

  const Cycle multigrid_cycle0 = Cycle::Z();
  const Cycle multigrid_cycle1 = Cycle::F();

  /* vof */
  const CurvMethod curv_method = CurvMethod::HF();
  const CurvMethod wall_curv_method = CurvMethod::HFnormalXZ();

  const TopoMethod topo_method = TopoMethod::Hybrid();

  const bool detachment_model = false;
  const bool subgrid_method = true; /* use slic subgrid */
  const bool use_fs_interp = false;
  const bool store_pressure_extrap = false;
  const int niter_pressure_extrap = 1000;

  /* enthalpy */
  const bool use_ht_resistance = false;//true;

  /* phase change - 4 version */
  const bool use_second_order_accuracy = true;
  const bool discard_points_near_interface = false;//true;
  const bool use_unconditional_extrapolation = false;

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();
