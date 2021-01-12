/******************************************************************************/
/* ------------ boundary or initial conditions */

  /* temperatures */
  const real tsat0 = 0.0;
  const real twall = tsat0 + deltat_wall;
  const real tout  = tsat0 + deltat_out;
  const real tnucl = tsat0 + deltat_nucl;

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
  real t_per_plot = 5e-6;
  if(case_flag==0)
    t_per_plot = 0.1;

  /* steps per backup */
  const int n_per_backup = 5000;

  /* simple alg */
  const int mSimple = 1;

  /* dt settings */
  const real initdtcoef = 1./50.;

  /* only liquid beyond this (fractional height) */
  const real zmax_mult = 0.9;

/* ------------ optional simulation settings */

  /* multigrid */
  const bool multigrid_stop_if_diverging = true;
  const bool multigrid_use_linf = false;//true;
  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 20;

  MaxIter multigrid_mm_smooth1 = MaxIter(35);
  MaxIter multigrid_mm_smooth2 = MaxIter(40);
  MaxIter multigrid_mm_solve = MaxIter(110);
  MaxIter multigrid_mm_stale1 = MaxIter(15);
  MaxIter multigrid_mm_stale2 = MaxIter(-1);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm_smooth1,multigrid_mm_smooth2,multigrid_mm_solve};
  std::array<MaxIter,3> multigrid_mstale = {multigrid_mm_stale1,multigrid_mm_stale1,multigrid_mm_stale2};

  ResRat multigrid_rr = ResRat(-1.);
  ResTol multigrid_rt = ResTol(3e-5);

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
  const ConvScheme cs_enth = ConvScheme::superbee();
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
