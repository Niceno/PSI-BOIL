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
  const real surftens_dt_coef = 1.;
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
  const bool multigrid_use_linf = false;//true;
  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 20;

  MaxIter multigrid_mm_smooth = MaxIter(20);
  MaxIter multigrid_mm_solve = MaxIter(100);
  MaxIter multigrid_mm_stale = MaxIter(-1);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm_smooth,multigrid_mm_smooth,multigrid_mm_solve};
  std::array<MaxIter,3> multigrid_mstale = {multigrid_mm_stale,multigrid_mm_stale,multigrid_mm_stale};

  ResRat multigrid_rr = ResRat(-1.);
  ResTol multigrid_rt = ResTol(5e-5);

  const Cycle multigrid_cycle0 = Cycle::Z();
  const Cycle multigrid_cycle1 = Cycle::F();

  /* vof */
  const CurvMethod curv_method = CurvMethod::HF();
  const CurvMethod wall_curv_method = CurvMethod::HFnormalXZ();

  const AdvectionMethod advect_method = AdvectionMethod::BoundedSplit();
  //const AdvectionMethod advect_method = AdvectionMethod::NaiveSplit();
  const TopoMethod topo_method = TopoMethod::Hybrid();

  const bool detachment_model = false;
  const bool subgrid_method = true; /* use slic subgrid */
  const bool use_fs_interp = false;
  const bool store_pressure_extrap = false;
  const int niter_pressure_extrap = 1000;

  /* enthalpy */
  const AccuracyOrder ao_efd_conv = AccuracyOrder::First();//Second();//Third();
  //const AccuracyOrder ao_efd_conv = AccuracyOrder::Third();
  const bool use_wall_resistance = true;

  /* phase change - 4 version */
  const AccuracyOrder ao_pc = AccuracyOrder::FourthUpwind();
  const bool discard_points_near_interface = false;//true;
  const bool use_unconditional_extrapolation = false;

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();
