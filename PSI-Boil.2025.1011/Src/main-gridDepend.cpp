#include "Include/psi-boil.h"
#include <iomanip>
#include <string>
#include <cstring>
#include "update_step.cpp"
using namespace std;

const int  gLevel = 32; //8

/* Domain */
const int  NX1 = 16*gLevel/8;
const int  NX2 =  2*gLevel/8;
const int  NY  = 8*gLevel/8;
const int  NZ  = 12*gLevel/16;  // 16 -> 24
const int  NmZ =  2*gLevel/8;  // 8  -> 
const real LY  =  0.005/8.0;
const real LX1 =  0.01/8.0;
const real LX3 =  0.012/8.0;
const real LZ  =  0.00036;
const real LmZ = -0.001/32.0;

/* parameter for boundary conditions */
const real prs = 10.5*1e+5;  // system pressure (Pa)
const real tsat0 = 0.0;
const real dtsub = 10.0;
const real qflux = 1.0e+6;  // heater power (W/m2)
const real thickITO = 0.7e-6; // thickness of ITO heater (m)
const real qsrc = qflux/thickITO;  // heater power (W/m3)

/* constants */
const real gravity = 9.8;

/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc==1){
    boil::oout<<"One command line argument is required!"<<"\n";
    boil::oout<<"./Boil wmin (wmin::wall time in minute)"<<"\n";
    exit(0);
  }
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  const real dx = LX1/real(NX1);
  Grid1D gx1( Range<real>(0.0, LX1), NX1, Periodic::no());
  Grid1D gx2( Range<real>(LX1, LX3), Range<real>(1.2*dx,2*dx), NX2, Periodic::no());
  Grid1D gx( gx1, gx2, Periodic::no());
  Grid1D gy( Range<real>(-LY/2.0, LY/2.0), NY, Periodic::yes());
  boil::oout<<"Grid1D: gy\n";

  Grid1D gz1( Range<real>(LmZ,0.0), Range<real>(-LmZ/NmZ*2.0,thickITO), NmZ, Periodic::no());
  boil::oout<<"Grid1D: gz1\n";
  Grid1D gz2( Range<real>(0.0,LZ), Range<real>(dx/4.0,dx*1.2),NZ, Periodic::no());
  //Grid1D gz2( Range<real>(0.0,0.5*LZ), Range<real>(dx/4.0,dx),NZ*3/4, Periodic::no());
  //boil::oout<<"Grid1D: gz2 NZ="<<NZ<<"\n";
  //Grid1D gz3( Range<real>(0.5*LZ,LZ),  Range<real>(1.2*dx,4.0*dx), NZ/4, Periodic::no());
  //boil::oout<<"Grid1D: gz3\n";
  //Grid1D gztmp( gz1, gz2, Periodic::no());
  //Grid1D gz( gztmp, gz3, Periodic::no(), BndGrid::wall(), BndGrid::symmetry());
  Grid1D gz( gz1, gz2, Periodic::no(), BndGrid::wall(), BndGrid::symmetry());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Domain dom(gx, gy, gz, & floor);
  int nig=gx.ncell();
  int njg=gy.ncell();
  int nkg=gz.ncell();
  //boil::plot->plot(dom,"dom");
  dom.save("grid32.grd");

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(dom), xyz(dom);             // velocity
  Scalar press(dom), p  (dom), f  (dom); // pressure
  Scalar tpr(dom), q  (dom);             // temperature
  Scalar c  (dom), g  (dom), step(dom), kappa(dom); // color function
  Scalar mdot  (dom); // phase change mdot 2
  Scalar wd(dom), mu_t(dom), yplus(dom);// eddy viscosity
  Scalar sdummy(dom);  // temporary or nousage
  Scalar id_bubble(dom);  // for Floodfill

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    //uvw.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    //uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet(), 0.0, 0.0, 0.5 ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell()),
              Range<int>(1,NmZ), BndType::inlet(), 0.0, 0.0, 0.0 ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell()), Range<int>(NmZ+1,gz.ncell()), 
              BndType::inlet(), "1.3201*(1.0-((0.00589-z)/0.00589)^1.2)^0.2", 0.0, 0.0 ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  }
  press.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );

  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  f    = p.shape();
  q    = p.shape();
  mdot = p.shape();
  kappa = p.shape();
  g     = p.shape();
  yplus = p.shape();
  mu_t  = p.shape();
  wd    = p.shape();
  sdummy  = p.shape();
  id_bubble = p.shape();

  tpr.bc().add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell())
            , Range<int>(1,NmZ), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::imin(), Range<int>(1, gy.ncell())
            , Range<int>(NmZ+1,gz.ncell()), BndType::dirichlet(), tsat0-dtsub) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );

  c.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(), 1.0 ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
  step = c.shape();

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(dom), liquid(dom), sapphire(dom); // zircaloy(dom)
  //zircaloy.lambda  (13.94);
  //zircaloy.rho    (6537.);
  //zircaloy.cp     (295.9*6537.);
  sapphire.rho    (3980.0);
  sapphire.cp     (750.0*3980.0);
  sapphire.lambda (35.0);

  const real tsat = IF97::Tsat97(prs);
  vapor.mu(IF97::viscvap_p(prs));
  vapor.rho(IF97::rhovap_p(prs));
  vapor.cp(IF97::cpvap_p(prs)*IF97::rhovap_p(prs));
  vapor.lambda(IF97::tcondvap_p(prs));

  liquid.mu(IF97::viscliq_p(prs));
  liquid.rho(IF97::rholiq_p(prs));
  liquid.cp(IF97::cpliq_p(prs)*IF97::rholiq_p(prs));
  liquid.lambda(IF97::tcondliq_p(prs));
  Matter mixed(liquid, vapor, &step);
  mixed.sigma(IF97::sigma97(tsat));
  const real latent=IF97::hvap_p(prs)-IF97::hliq_p(prs);

  /* estimate betal */
  const real liquid_drhodt = (IF97::rholiq_p(1.01*prs)-IF97::rholiq_p(0.99*prs))
                             /(IF97::Tsat97(1.01*prs)-IF97::Tsat97(0.99*prs));
  const real vapor_drhodt=0.0; //[kg/m3K]

  boil::oout<<"properties at pressure "<<prs<<boil::endl;
  boil::oout<<"Tsat (K)= "<<tsat<<" vapor= "<<IF97::viscvap_p(prs)<<" (Pa.s) "
           <<IF97::rhovap_p(prs)<<" (kg/m3) "<<IF97::cpvap_p(prs)<<" (J/kg.K) "
	   <<IF97::tcondvap_p(prs)<<" (W/m.K)\n";
  boil::oout<<"Tsat (C)= "<<tsat-273.15<<" liquid= "<<IF97::viscliq_p(prs)<<" (Pa.s) "
           <<IF97::rholiq_p(prs)<<" (kg/m3) "<<IF97::cpliq_p(prs)<<" (J/kg.K) "
           <<IF97::tcondliq_p(prs)<<" (W/m.K)\n";
  boil::oout<<"sigma= "<<IF97::sigma97(tsat)<<" latent= "<<latent
           <<" drho/dT= "<<liquid_drhodt<<"\n";

  const real rhol  = liquid.rho()->value();
  const real rhov  = vapor.rho()->value();
  const real mul   = liquid.mu()->value();
  const real muv   = vapor.mu()->value();
  const real sigma = mixed.sigma()->value();
  real const mu_t_max =100.0*std::max(mul,muv);  // previously 10.0

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int  ndt = 50000000;
  const real tint = 1.0e-4;  // time interval for dat-output
  const int  nint = 10000;   // time interval for bck-output
  const real dxmin = dom.dxyz_min(); // minimum size of cell in liquid
#if 1
  const real dt  = 10.0*pow(vapor.rho()->value()*pow(dxmin,3.0)
	             / (2.0*3.1415*sigma),0.5);
#else
  const real dt = 8e-7;
#endif
  const real cfl_limit=0.25;
  Times time(ndt, dt);
  time.print_time(false);
  time.set_coef_dec(0.5);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(dom, Prec::ic2());
  //Krylov * solver2 = new BiCGS(dom, Prec::di());
  Krylov * solver2 = new CG(dom, Prec::di());  // Change on 2024.1128 because 
					       // of floating point exception

  Pressure pr( p,   f,   uvw, time, solver, &mixed );
  Momentum ns( uvw, xyz,        time, solver, &mixed );
  ns.convection_set(TimeScheme::forward_euler());
  EnthalpyFD en(tpr, q, c, uvw, time, solver2, & mixed, tsat0, & sapphire);
  en.convection_set(TimeScheme::forward_euler());
  en.diffusion_set(TimeScheme::backward_euler());
  CIPCSL2 conc (c,  g, kappa, uvw, time, solver);
  conc.set_itsharpen(10);
  conc.set_globalSharpen();
  conc.set_nredist(1);
  conc.set_cangle(0.0);

  const real rseed = 1.2 * dx;
  const real dmicro_min = 1.0e-10; // nonuse
  Nucleation nucl( &c, &tpr, &q, &time, sdummy, &mixed,
    rseed, dmicro_min, latent, conc.get_cangle(), tsat0, &sapphire);
  nucl.set_seed_period(1e-5);    // continue to plant color function
  nucl.set_pre_heat_sink(true);  // give heat sink before plant of color function
  nucl.set_micro_exists(false);  // micro-layer model is not used
  nucl.set_rseed_plus(2.0);      // clr detect range for replant

  PhaseChange pc(mdot, tpr, q, c, g, f, step, uvw,
                 time, &mixed, latent, tsat0, &sapphire, &nucl);

  AC sol( &pr );
  sol.stop_if_diverging(true);
  sol.min_cycles(3);
  sol.max_cycles(10);

  /* wall distance */
  Distance di(wd, sdummy, uvw, time, solver);

  /* eddy viscosity */
  Model tm;

  /* track bubble ID */
  Floodfill flood(c, id_bubble, &uvw, time);
  flood.set_out_freq(10);
  boil::oout<<"main:Floodfill:output frequency= "<<flood.get_out_freq()<<"\n";
  flood.set_smallest_cvol(1);
  flood.set_range(Range<real>(0.0,LX1),Range<real>(-0.5*LY,0.5*LY),Range<real>(0.0,LZ));

  /*-------------------+
  |  check if restart  |
  +-------------------*/

  std::fstream input;
  int irun = 0;
  if(boil::cart.iam()==0){
    input.open("run.txt", std::ios::in);
    if( !input.fail() ) {
      input >> irun;
      boil::oout<<"read irun.  irun= "<<irun<<"\n";
    }
    input.close();
  }
  boil::cart.sum_int(&irun);
  if (irun==1){
    boil::oout<<"exit job due to irun=1"<<"\n";
    exit(0);
  }

  if(boil::cart.iam()==0){
    std::fstream output;
    output.open("run.txt", std::ios::out);
    output << 1 << boil::endl;
    output.close();
  }

  int ts=0;
  bool restart = false;

  input.open("time.txt", std::ios::in);
  if( !input.fail() ) {
    restart=true;
  }

  if( restart ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
  }
 
  if( restart ) {
    uvw   .load("uvw", ts);
    press .load("press", ts);
    tpr   .load("tpr", ts);
#if 1
    /* normal start */
    conc  .load("conc", ts);
    nucl  .load("nucl", ts);
    wd    .load("wd",   ts);
#if 0
    flood.load("flood", ts);
#endif
#else
    /* after changeGrid */
    /* load conc-phi as a scalar */ 
    c.load("conc-phi", ts);
    c.exchange_all();
    /* calculate conc-sigma etc. */ 
    conc.init();
    /* zsite can be changed */
    real zsite=rseed*cos(conc.get_cangle()/180.0*boil::pi);
    boil::oout<<"main:change_site_z= "<<zsite<<" "<<rseed<<"\n";
    nucl  .load("nucl", ts, &zsite);
    /* distance function from wall */
    wd    .load("wd",   ts);  // load as an initial value
    sdummy = 0.0;             // reset source term
    di.compute(1.0e-15);
    boil::plot->plot(wd,"wd",ts);
#endif
    /* calculate step in case of restart */
    update_step(c, step, sdummy);
  }

  else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    /* compute distance function */
    sdummy = 0.0; // reset source term for di
    di.compute(1.0e-15);
    //boil::plot->plot(wd,"wd",0);

    /* velocity */
    Comp m = Comp::u();
    for_vmijk(uvw,m,i,j,k){
      real zz = uvw.zc(m,k);
        if(zz<0.0) {
          uvw[m][i][j][k] = 0.0;
        } else {
          uvw[m][i][j][k] = 1.3201*(pow(1.0-pow(((0.00589-zz)/0.00589),1.2),0.2));
        }
    }

    /* color function */
    for_vijk(c,i,j,k)
      c[i][j][k] = 1.0;
    c.exchange_all();
    update_step(c, step, sdummy);
    conc.init();

    /* temperature */
    tpr = tsat0-dtsub;
    for_vijk(tpr,i,j,k){
      if(tpr.zc(k)<0.0){
        tpr[i][j][k] = tsat0;
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

    /* set nucleation sites */
#if 1
    ifstream fin("site-G1000.txt");
    if( fin.fail() ) {
      boil::oout<<"Cannot open nucleation site file\n";
      exit(0); 
    } else {
      boil::oout<<"Open nucleation site file\n";
      string ss; 
      getline(fin, ss); // skip first line
      real nsx,nsy,Tact;
      real zsite=rseed*cos(conc.get_cangle()/180.0*boil::pi);
      boil::oout<<"main:zsite= "<<zsite<<" rseed= "<<rseed<<"\n";
      int ins=0;
      std::string line;
      while (std::getline(fin, line)){
        std::istringstream iss(line);
        real nsx, nsy, Tact;
        iss >> nsx >> nsy >> Tact;
        if(0.02*LX1<nsx && nsx<LX1 && -0.5*LY<nsy && nsy<0.5*LY) {
          nucl.add(Site(nsx, nsy, zsite, Tact, -0.002));
          if(ins%100==0)
            boil::oout<<"main:nucl.add= "<<ins<<" "<<nsx<<" "<<nsy<<" "<<Tact<<"\n";
          ins++;
	}
      	// -0.002 is nonuse: bottom of previous bubble
        //if(ins%1000==0){
        //  boil::oout<<"main:nucl.add= "<<ins<<" "<<nsx<<" "<<nsy<<" "<<Tact<<"\n";
        //}
      }
    }
    fin.close();
#else
    #include "nsd.cpp"
#endif
    /* plot initial condition */
    boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,yplus,id_bubble,
                    "uvw-c-tpr-press-mdot-mut-yplus-id",0);
  }
  input.close();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    boil::oout << "########################" << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "# DT:        " << time.dt() << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() 
               << "/"             << time.total_steps() << boil::endl;
    boil::oout << "# WTIME:     " << boil::timer.current_min() << boil::endl;
    boil::oout << "########################" << boil::endl;

    /*---------------------+
    |  reset source terms  |
    +---------------------*/
    /* body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /* heat source */
    for_vijk(q,i,j,k)
      q[i][j][k] = 0.0;

    for_vk(q,k){
      if(-thickITO<q.zc(k) && q.zc(k)<0.0){
        for_vij(q,i,j){
          q[i][j][k] = qsrc*q.dV(i,j,k);  // [W/m3]*[m3]=[W]
        }
      }
    }

    /*------------------+
    |  turbulence model |
    +------------------*/
    /* Smagorinsky model */
    //tm.smagorinsky( & ns, & mu_t, 0.17);  // without damping
    const real dist_max = 1.0e-3; // at z = 1mm, yplus ~ 180. damping ~ 0 if yplus > 40.
    tm.smagorinsky( & ns, & mu_t, 0.17, &dist_max, &wd, &yplus);  // with damping
    /* WALE model */
    //tm.wale( & ns, & mu_t, 0.1);  // 0.325^2 ~ 0.544^2
    /* viscous term in wall adjacent cells */
    tm.tau_wall( & ns, wd, & xyz );

    /* limit eddy viscosity */
    for_avijk(mu_t,i,j,k) {
      mu_t[i][j][k] = std::min(mu_t[i][j][k],mu_t_max);
    }

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update( &mu_t );
    //pc.micro(&xyz);
    update_step(c, step, sdummy);  // 0119 need to update step for with IB
    ns.vol_phase_change(&f);

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* gravity force */
    Comp m = Comp::u();
    for_vmijk(xyz,m,i,j,k){
      real phil=step[i][j][k];
      real phiv=1.0-phil;
      real deltmp = tpr[i][j][k]-tsat0;
      real rhomix = (rhol + liquid_drhodt*deltmp)*phil
                  + (rhov + vapor_drhodt*deltmp)*phiv;
      if(dom.ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }
    xyz.exchange();

    /* surface tension */
    conc.tension(&xyz, mixed, step);

#if 1
    /* increase viscosity in outlet region */  // COPY
    const real x0 = LX1;
    const real x1 = LX3;
    for_avi(c,i){
      if(c.xc(i)>x0){
        real coef=std::min((c.xc(i)-x0)/(x1-x0),1.0);
        for_avjk(c,j,k){
          mu_t[i][j][k] = coef * liquid.mu()->value() * 10;
        }
      }
    }
#endif
     
    /* intermediate velocity */
    ns.discretize( &mu_t );

    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1.0e-32));

    /* pressure poisson equation */
    p = 0.0;
    sol.vcycle(ResRat(1.0e-16));

    /* update velocity and pressure */
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        if(pmin>press[i][j][k]) pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k){
      if(dom.ibody().on(i,j,k)){
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }
    press.bnd_update();
    press.exchange_all();

#if 1
    /* limit velocity */
    for_m(m)
      for_avmijk(uvw,m,i,j,k) {
        uvw[m][i][j][k] = max(-25.0,min(25.0,uvw[m][i][j][k]));
      }
#endif

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc.advance();
    conc.totalvol(Range<real>(0.0, LX1),
                  Range<real>(-LY/2.0, LY/2.0),
                  Range<real>(0.0, LZ));

    /*---------------------------+
    |  replant seed or cut neck  |
    +---------------------------*/ 
    for(int nsd=0; nsd<nucl.size(); nsd++){
      nucl.sites[nsd].set_allow_replant(true);
    }
    nucl.replant();

    /*------------------------------+
    |  outlet region: delete liquid |
    +------------------------------*/
    const real x0b = LX1 + 0.5*(LX3-LX1);
    const real x1b = LX3;
    for_avi(c,i){
      if(c.xc(i)>x0b){
        real coef=std::min((c.xc(i)-x0b)/(x1b-x0b),1.0);
        for_avjk(c,j,k){
          c[i][j][k]= (1.0-coef)*c[i][j][k] + coef*0.0;
        }
      }
    }
    /* update clr in cipcsl2 after seed, cutneck and outlet-region */
    c.bnd_update();
    c.exchange_all();
    conc.update_node(c);
    update_step(c, step, sdummy);

    /* output min & max of color function */
    conc.color_minmax();
    boil::oout<<"main:color_min,max= "<<time.current_time()<<" "
             <<conc.minval()<<" "<<conc.maxval()<<"\n";

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    en.discretize( &mu_t );
    en.new_time_step( &mu_t );
    en.solve(ResRat(1e-16),"enth");
    //enthFD.solve_sor(24, 1.2, "enthFD");

    /* outlet region */
    for_vi(tpr,i){
      if(c.xc(i)>LX1+(LX3-LX1)*0.5){
          for_vk(tpr,k) {
            if (tpr.zc(k)<0) {
              for_vj(tpr,j) {
                tpr[i][j][k]=tpr[i-1][j][k];
              }
            } else {
              for_vj(tpr,j) {
                tpr[i][j][k]=max(tpr[i-1][j][k],tsat0);
              }
            }
         }
      }
    }
    tpr.bnd_update();
    tpr.exchange_all();

    /*------------+
    |  floodfill  |
    +------------*/
    boil::oout<<"main:floodfill\n";
    flood.identify_regions();

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax = ns.cfl_max();
    boil::oout<<"main:cflmax= "<<time.current_time()<<" "<<cflmax<<"\n";
    /* interface is always included because of outlet */
    time.control_dt(cflmax, cfl_limit, dt);

    if (time.dt()<1e-8) {
      boil::oout<<"Too small time step: "<<time.dt()<<"\n";
      exit(0);
    }

    /*--------------+
    |  output data  |
    +--------------*/
    /* heat flux and area */
    en.hflux_wall_ib( Range<real>(0.0, LX1), Range<real>(-LY/2.0, LY/2.0),&mu_t);

    /* average temperature in ITO */
    real sum_tcb=0.0;
    real vol_tcb=0.0;
    int n_tc=0;
    for_vk(c,k) {
      if(-thickITO<c.zc(k) && c.zc(k)<0.0){
        for_vij(c,i,j) {
          if(c.xc(i)<LX1) {
            sum_tcb+=tpr[i][j][k]*tpr.dV(i,j,k);
            vol_tcb+=tpr.dV(i,j,k);
          }
        }
      }
    }
    boil::cart.sum_real(&sum_tcb);
    boil::cart.sum_real(&vol_tcb);

    std::cout.setf(std::ios_base::scientific);
    std::cout<< std::setprecision(12);
    boil::oout<<"twall= "<<time.current_time()<<" "<<sum_tcb/vol_tcb<<"\n";
    std::cout<< std::setprecision(8);
    std::cout.unsetf(std::ios_base::floatfield);

    /* tecplot files */
    if((time.current_time()) / (tint) >= real(iint) ) {
      iint = int(time.current_time() / tint);
      tpr.exchange_all();
      boil::plot->plot(uvw,c,tpr,press,mdot,mu_t,yplus,id_bubble,
                      "uvw-c-tpr-press-mdot-mut-yplus-id",iint);
      iint = int(time.current_time()/tint) + 1;
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % nint==0) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
      nucl .save("nucl",   time.current_step());
      wd   .save("wd",   time.current_step());
      flood.save("flood",time.current_step());
    }


    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      conc .save("conc",  time.current_step());
      tpr  .save("tpr",   time.current_step());
      nucl .save("nucl",   time.current_step());
      wd   .save("wd",   time.current_step());
      flood.save("flood",time.current_step());

      uvw  .rm("uvw",   ts);
      press.rm("press", ts);
      conc .rm("conc",  ts);
      tpr  .rm("tpr",   ts);
      nucl .rm("nucl", ts);
      wd   .rm("wd", ts);
      flood.rm("flood", ts);
    }

    
    if((time.current_step()) % (nint)==0 ) {
      if( boil::cart.iam()==0) {
        std::fstream output;
        std::stringstream ss;
        ss <<"time-"<<time.current_step()<<".txt";
        std::string fname = ss.str();
        int len = fname.length();
        char * cfname = new char[len+1];
        memcpy(cfname, fname.c_str(), len+1);
        output << std::setprecision(16);
        output.open(cfname, std::ios::out);
        output << time.current_step() << boil::endl;
        output << time.current_time()+time.dt() << boil::endl;
        output << time.dt() << boil::endl;
        output.close();
      }
    }
 
    if( boil::timer.current_min() > (wmin)
      || time.current_step()==time.total_steps()) {
      std::fstream output;
      output << std::setprecision(16);
      output.open("time.txt", std::ios::out);
      output << time.current_step() << boil::endl;
      output << time.current_time()+time.dt() << boil::endl;
      output << time.dt() << boil::endl;
      output.close();
      output.open("run.txt", std::ios::out);
      output << 0 << boil::endl;
      output.close();
      boil::timer.stop();
      boil::timer.report();
      exit(0); 
    }
  }

  boil::oout << "finished" << boil::endl;
  boil::timer.stop();
  boil::timer.report();

}
/*-----------------------------------------------------------------------------+
 '$Id: main-dam.cpp,v 1.12 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
