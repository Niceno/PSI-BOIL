//#define USE_PHASE_VEL_EFD
#define USE_SOLID
//#define SETUP_ONLY
#define USE_BIG
#define USE_PC4
#define CASE 0
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

  if(argc<7){
    boil::oout<<"Six command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin glevel gstage cangle deltat pressure[bar]"<<"\n";

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

  const real prs = 101325.*atof(argv[6]); /* pressure */
  boil::oout<<"pressure= "<<prs<<"\n";

/******************************************************************************/
/* ------------ rescaling factors */
  const real xmult = 1e0;
  const real tmult = 1e0;
  const real mmult = xmult*tmult;

/******************************************************************************/
/* ------------ boundary or initial conditions */
  const real tsat0 = 0.0;

  const real twall = tsat0 + deltat;
  const real tsat0_K = IF97::Tsat97(prs);

  const real twall0 = twall;
  const real tseed = twall0-0.001;

#if CASE == 0
  #if 0
  const real tout = tsat0;
  #else
  const real tout = twall0;//tsat0;
  #endif
#else
  const real tout = tsat0;
#endif

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
  const real Mv = IF97::get_MW();
  const real muv = IF97::viscvap_p(prs);
  const real rhov = IF97::rhovap_p(prs);
  const real cpv = IF97::cpvap_p(prs)*rhov;
  const real lambdav = IF97::tcondvap_p(prs);

  const real mul = IF97::viscliq_p(prs);
  const real rhol = IF97::rholiq_p(prs);
  const real cpl = IF97::cpliq_p(prs)*rhol;
  const real lambdal = IF97::tcondliq_p(prs);

  const real sig = IF97::sigma97(tsat0_K);
  const real latent=IF97::hvap_p(prs)-IF97::hliq_p(prs);

  const real betal = (7.03+(tsat0_K-273.15-100.)/160.*(22.1-7.03))*1e-4; /* roughly */
  const real betav = 1./tsat0_K; /* ideal gas approximation */

  const real alpl = lambdal/cpl;
  const real alpv = lambdav/cpv;

  boil::oout<<"properties at pressure "<<prs<<boil::endl;
  boil::oout<<"vapprop= "<<Mv<<" "<<muv<<" "<<rhov<<" "
                         <<cpv/rhov<<" "<<lambdav<<" "<<betav<<" "<<alpv<<boil::endl;
  boil::oout<<"liqprop= "<<Mv<<" "<<mul<<" "<<rhol<<" "
                         <<cpl/rhol<<" "<<lambdal<<" "<<betal<<" "<<alpl<<boil::endl;
  boil::oout<<"twoprop= "<<tsat0_K<<" "<<sig<<" "<<latent<<boil::endl;

  const real Jal = cpl*(twall-tsat0)/(latent*rhov);
  boil::oout << "Jal= "<<Jal<<boil::endl;

  /* heater */
  /* sapphire */
  const real rhosol = 3980.0;
  const real cpsol = 750*rhosol;
  const real lambdasol = 35.;

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

  const real correl_reduc = 9.*boil::pi/8.*mul*alpl/sig*Jal*Jal/(R0mult*DX0);
  boil::oout<<"correl: theta^3/lnS= "<<correl_reduc<<" rad^3"<<boil::endl;

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
