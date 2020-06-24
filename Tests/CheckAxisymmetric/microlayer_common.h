//#define USE_PHASE_VEL_EFD
#define USE_SOLID
//#define SETUP_ONLY
#define USE_BIG
#define USE_PC4
#define CASE 1
/*  0: Toulouse-like
    1: Duan [MIT] (2013)
 */

#if CASE == 0
  #define USE_BOTTOM_DIRICHLET
#endif
#ifndef USE_SOLID
  #define USE_BOTTOM_DIRICHLET
#endif

/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc<6){
    boil::oout<<"Five command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin glevel gstage cangle deltat"<<"\n";

    exit(0);
  }

/******************************************************************************/
/* ------------ input from command line */
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const int gLevel = atoi(argv[2]); /* domain dimensions */
  boil::oout<<"glevel= "<<gLevel<<"\n";

  const int gStage = atoi(argv[3]); /* domain dimensions */
  boil::oout<<"gstage= "<<gStage<<"\n";

  const real cangle = atof(argv[4]); /* contact angle */
  boil::oout<<"cangle= "<<cangle<<"\n";

  const real deltat = atof(argv[5]); /* superheat */
  boil::oout<<"deltat= "<<deltat<<"\n";

/******************************************************************************/
/* ------------ rescaling factors */
  const real xmult = 1e0;
  const real tmult = 1e0;
  const real mmult = xmult*tmult;

/******************************************************************************/
/* ------------ boundary or initial conditions */
  const real tsat0 = 0.0;
  const real tout = tsat0;

  const real twall = tsat0 + deltat;
  const real tsat0_K = 373.12;

  const real twall0 = twall;
  const real tseed = twall0-0.001;

/******************************************************************************/
/* ------------ values to be directly rescaled */
  real gravity = boil::g; /* m/s2 */
  real R = boil::R;
  
#if CASE == 1
  real qflux=28.7e3;
#else
  real qflux=0.;
#endif
  real seedper = 1e-5; /* s */

  /* rescaling */
  gravity *= xmult/tmult/tmult; /* [m/s2] */
  R *= mmult*xmult*xmult/tmult/tmult; /* [kgm2/s2/mol/K] */

  seedper *= tmult;

  qflux *= mmult/tmult/tmult/tmult;     // [kg/s3]

/******************************************************************************/
/* ------------ numerical simulation settings */

  /* total number of steps */
  const int ndt = 10e6; /* inconsequential */

  /* total time */
  const real tend = 1e-3 * tmult;

  /* steps per backup */
  const int n_per_backup = 5000;

  /* if yes, plotting each t_per_plot seconds. Else, each n_per_plot steps */
  const bool use_t_per_plot = true;
  const real t_per_plot = tend/100.;
  const int n_per_plot = ndt/200;

  /* dt settings */
  const real surftens_dt_coef = 10.;
  const real initdtcoef = 1./10.;

  /* cfl with and without interfaces */
  const real cfl_with = 0.1; //0.05
  const real cfl_wo   = 0.2;
 
  /* radius of init bubble in terms of dxmin */
  const real R0mult = 20.;

/* ------------ optional simulation settings */

  /* multigrid */
  const bool multigrid_stop_if_diverging = true;
  //const bool multigrid_stop_if_diverging = false;

  const int multigrid_min_cycles = 1;
  const int multigrid_max_cycles = 10+2*gLevel;

  const int multigrid_niter = 30;
  MaxIter multigrid_mm = MaxIter(multigrid_niter);
  std::array<MaxIter,3> multigrid_mi = {multigrid_mm,multigrid_mm,multigrid_mm};

  ResRat multigrid_rr = ResRat(5e-5);

  const Cycle multigrid_cycle0 = Cycle::Z();
  const Cycle multigrid_cycle1 = Cycle::F();


  /* vof */
  const CurvMethod curv_method = CurvMethod::HF();
  const CurvMethod wall_curv_method = CurvMethod::HFnormalXZ();
  //const CurvMethod wall_curv_method = CurvMethod::HFmixedXZ();
  const int Nfilm_crit = 4;
  const TopoMethod topo_method = TopoMethod::Hybrid();

  const bool detachment_model = false;
  const bool subgrid_method = true; /* use slic subgrid */
  const bool use_fs_interp = false;
  const bool store_pressure_extrap = false;
  const int niter_pressure_extrap = 1000;

  /* phase change - VOF version */
  const bool near_wall_modelling = true; /* use near wall modelling */
  /* phase change - 4 version */
  const bool use_second_order_accuracy = true;
  const bool discard_points_near_interface = true;
  const bool use_unconditional_extrapolation = false;

/******************************************************************************/
/* ------------ material properties */
  const real Mv = 18.015e-3;
  const real muv = 1.228e-5;
  const real rhov = 0.5974;
  const real cpv = 2034*rhov;
  const real lambdav = 0.024;

  const real mul = 2.82e-4;
  const real rhol = 958;
  const real cpl = 4216*rhol;
  const real lambdal = 0.677;

  const real sig = 0.058;
  const real latent=2256e3;

  const real betal = 7.52e-4;
  const real betav = 1./tsat0_K; /* ideal gas approximation */

  /* heater */
  /* sapphire */
  const real rhosol = 3980.0;
  const real cpsol = 750*rhosol;
  const real lambdasol = 35.;

  const real Jal = cpl*(twall-tsat0)/(latent*rhov);
  boil::oout << "Jal= "<<Jal<<boil::endl;

/******************************************************************************/
/* ------------ domain dimensions */
  real L0 = 0.25e-3 * gStage;
  const int N0 = 512*gLevel;

  const int ARx = 2;
  const int NX1 = ARx*N0;
  real LX1 = real(ARx)*L0;
#if CASE == 0
  const int ARz = 3;
  /* in experiment, ITO: 700 nm, CaF2: 10 mm;
     here, CaF2: 4*DX0 */
  const int NZ0 = 4;
#elif CASE == 1
  const int ARz = 3;
#endif

  LX1 *= xmult;

/******************************************************************************/
/* ------------ calculated values */
  const real DX0 = LX1/real(NX1);
  boil::oout<<"DX= "<<DX0<<boil::endl;
  const real DZ0 = DX0;

#if CASE == 1
  /* in experiment, ITO: 700 nm, sapphire: 250 um;
     here, sapphire: ~100 um; note: proc and solid boundary mustnt overlap! */
  //const int NZ0 = N0/(gStage*5);
  const int NZ0 = N0/gStage-3*gLevel;
  /* should be 0.7 um */
  const int NZheat = std::ceil(0.7e-6/DZ0);
  const real LZ0 = -DZ0*NZ0;
  const real LZheat = -NZheat*DZ0;
  boil::oout<<"NZheat= "<<NZ0<<" "<<NZheat<<" "<<LZheat<<boil::endl;
#else
  const real LZ0 = -DZ0*NZ0;
  const real LZheat = LZ0;
#endif

  int NZ(0);
#ifdef USE_SOLID
  const int NZ1 = ARz*N0 - NZ0;
  NZ += NZ0+NZ1;
#else
  const int NZ1 = ARz*N0;
  NZ += NZ1;
#endif

  const real LZ1 = DZ0*NZ1;

  /* other parameters */
  const real zplant= LZ1/18.; // when bottom of bubble reaches zplant, next seed is set
  const real zmax=0.9*LZ1; /* only liquid beyond 0.9*LZ1 */

  /* heater power */
  real qsrc=qflux/fabs(LZheat);  // [W/m3]
  boil::oout<<"#qsrc= "<<qsrc<<" qflux= "<<qflux<<"\n";

/******************************************************************************/
/* ------------- setup finished */
/******************************************************************************/
/* below this line, NO VALUES AND SETTINGS CAN BE ENTERED! */
/******************************************************************************/
