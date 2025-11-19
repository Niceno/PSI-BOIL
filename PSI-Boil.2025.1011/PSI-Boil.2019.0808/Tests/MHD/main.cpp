/* Validation of Wang's experiment: rising argon bubble in GaInSn. 
  B = 1.969 (T), Eo = 2.44
  Z.H. Wang, S.D. Wang, X. Meng, M.J. Ni, UDV measurements of single bubble
  rising in a liquid metal Galinstan with a transverse magnetic field, 
  International Journal of Multiphase Flow, 94 (2017) 201-208. */
#include "Include/psi-boil.h"
#include <vector>
#define FLOODFILL

/* parameters */
const real radius = 0.004568/2.0 ;//Eo=2.44
const real dia = radius *2.0;

const real LX = 6.0*dia;
const real LY = LX;
const real LZ= 24*dia;

const int glevel = 6;
const int NX = 16*glevel;
const real dxmin = LX/real(NX); //to change
const int NY = NX;
const int NZ= LZ/dxmin;
const real z_init=2.0*dia;
const real gravity=9.8;
const real B0=1.969;

/****************************************************************************/
int main(int argc, char * argv[]) {

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

  /*--------+
  |  grids  |
  +--------*/
  Grid1D gx(Range<real>( -LX/2.0, LX/2.0 ), NX, Periodic::no());
  Grid1D gy(Range<real>( -LY/2.0, LY/2.0 ), NY, Periodic::no());
  Grid1D gz(Range<real>( 0.0,  LZ ), NZ, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);
	
  /*----------------+
  |  linear solver  |
  +----------------*/
  Krylov * solver = new BiCGS(d, Prec::di());
  Krylov * solverNS = new CG(d, Prec::di());

  /*------------------+
  |  define unknowns  |
  +------------------*/
   Vector uvw(d), xyz(d); // velocity
   Scalar press(d), p  (d), f  (d); // pressure
   Scalar c  (d), g  (d), kappa(d); // concentration
   Vector B(d),J(d) ; //magnetic flux density, electric current
   Scalar pot(d), pot_src(d) ; // electric potental
   Scalar id_bubble(d) ;

  /*-----------------------------+ 
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
}

  p.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press = p.shape();
  g=p.shape();
  kappa = p.shape();
  f=p.shape();
  pot=p.shape();
  pot_src=p.shape();
  id_bubble =p.shape();

  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );

  for_m(m) {
    B.bc(m).add( BndCnd( Dir::imin(), BndType::wall(), B0, 0.0, 0.0 ) );
    B.bc(m).add( BndCnd( Dir::imax(), BndType::wall(), B0, 0.0, 0.0 ) );
    B.bc(m).add( BndCnd( Dir::jmin(), BndType::wall(), B0, 0.0, 0.0 ) );
    B.bc(m).add( BndCnd( Dir::jmax(), BndType::wall(), B0, 0.0, 0.0 ) );
    B.bc(m).add( BndCnd( Dir::kmin(), BndType::wall(), B0, 0.0, 0.0 ) );
    B.bc(m).add( BndCnd( Dir::kmax(), BndType::wall(), B0, 0.0, 0.0 ) );
  }

  Comp m ;
  m = Comp::u();
    J.bc(m).add( BndCnd( Dir::imin(), BndType::wall(), 0.0 ,0.0, 0.0 ) );
    J.bc(m).add( BndCnd( Dir::imax(), BndType::wall(), 0.0, 0.0, 0.0 ) );
    J.bc(m).add( BndCnd( Dir::jmin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::jmax(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::kmin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::kmax(), BndType::neumann() ) );
    
  m = Comp::v(); 
    J.bc(m).add( BndCnd( Dir::imin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::imax(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::jmin(), BndType::wall(), 0.0, 0.0, 0.0 ) );
    J.bc(m).add( BndCnd( Dir::jmax(), BndType::wall(), 0.0, 0.0, 0.0 ) );
    J.bc(m).add( BndCnd( Dir::kmin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  m = Comp::w();
    J.bc(m).add( BndCnd( Dir::imin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::imax(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::jmin(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::jmax(), BndType::neumann() ) );
    J.bc(m).add( BndCnd( Dir::kmin(), BndType::wall(), 0.0, 0.0, 0.0 ) );
    J.bc(m).add( BndCnd( Dir::kmax(), BndType::wall(), 0.0, 0.0, 0.0 ) );
  
  /*----------------------+
  |  physical properties  |
  +----------------------*/
#if 0
  Matter nitrogen(d), mercury(d);
  nitrogen.mu      ( 1.77e-5 );
  nitrogen.rho     ( 1.17    );
  nitrogen.sigma_e ( 1.00e-15);
  mercury.mu       ( 1.50e-3 );
  mercury.rho      ( 1.35e+4 );
  mercury.sigma_e  ( 1.02e+6 );
  Matter mixed(mercury,nitrogen, &c);
  mixed.sigma(0.4535);
#endif
#if 1
  Matter argon(d), GaInSn(d);
  argon  .mu ( 1.176e-5 );
  argon  .rho( 1.654    );
  argon.sigma_e ( 1.00e-15);
  GaInSn.mu    ( 2.20e-3 );
  GaInSn.rho   ( 6.3615e3  );
  GaInSn.sigma_e  ( 3.27e6  );
  Matter mixed(GaInSn,argon, &c);
  mixed.sigma(0.533);
#endif

  /*-------+
  |  Time  |
  +-------*/
  //const real dt  = 5.0*pow(0.5*nitrogen.rho()->value()*pow(dxmin,3.0)/(2.0*3.1415*mixed.sigma()->value()),0.5);
  const real dt  = 5.0*pow(0.5*argon.rho()->value()*pow(dxmin,3.0)/(2.0*3.1415*mixed.sigma()->value()),0.5);
  boil::oout<<"dt= "<<dt<<"\n";
  const int nint=10000;
  Times time(200000, dt); // ndt, dt
  const real cfl_limit=0.25;
  const real tint=5.0e-2;

  /*-----------------+
  |  define equation  |
  +-----------------*/
  Pressure pr( p,   f,   uvw, time, solverNS, &mixed );
  AC multigrid( &pr );
  multigrid.stop_if_diverging(true);
  multigrid.min_cycles(3);
  multigrid.max_cycles(10);

  Momentum ns( uvw, xyz,      time, solverNS, &mixed );
  ns.diffusion_set (TimeScheme::backward_euler());
  ns.convection_set(TimeScheme::forward_euler());
  //ns.convection_set(ConvScheme::upwind()); 


  ElectPoten pt(pot, pot_src, B, J, uvw, time, solver, &mixed);
  AC multigrid_pt( &pt );
  multigrid_pt.stop_if_diverging(true);
  multigrid_pt.min_cycles(3);
  multigrid_pt.max_cycles(10);


  /*-----------------+
  | interface traking |
  +-----------------*/
  VOF conc(c, g, kappa, uvw, time, solverNS);
  conc.set_cangle(90.0);

  /*-------------------+
  | traking bubble ID  |
  +-------------------*/
#ifdef FLOODFILL
  Floodfill flood(c, id_bubble, &uvw, time);
  flood.set_out_freq(1);
  boil::oout<<"main:Floodfill:output frequency= "<<flood.get_out_freq()<<"\n";
#endif

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  std::fstream input;
  int irun = 0;
  if(boil::cart.iam()==0){
    input.open("run.txt", std::ios::in);
    if( !input.fail() ) {
      input >> irun;
      std::cout<<"read irun.  irun= "<<irun<<"\n";
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

 /*----------+
  |  restart  |
  +----------*/
  if( restart ) {
    real t,dtf;
    input >> ts;
    input >> t;
    input >> dtf;
    time.first_step(ts);
    time.current_time(t);
    time.set_dt(dtf);
    uvw   .load("uvw", ts);
    press .load("press", ts);
    c     .load("c", ts);
  } else {
    boil::oout << "######################" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "# START FROM SCRATCH #" << boil::endl;
    boil::oout << "#                    #" << boil::endl;
    boil::oout << "######################" << boil::endl;

    /*--------------------+
    |  initial condition  |
    +--------------------*/
    for_avijk(c,i,j,k) {
      c[i][j][k]=1.0;
    }

    // bubble creation 
    for_vijk(c,i,j,k) {
      real dist=sqrt(pow(c.xc(i),2.0)+pow(c.yc(j),2.0)+pow((c.zc(k)-z_init),2.0));
      if (dist<radius*0.75) {
        c[i][j][k]=0.0;
      } else if(dist<radius*1.25) {
        int mm=8;
        real x0=d.xn(i);
        real y0=d.yn(j);
        real z0=d.zn(k);
        real ddx=d.dxc(i)/real(mm);
        real ddy=d.dyc(j)/real(mm);
        real ddz=d.dzc(k)/real(mm);
        int itmp=0;
        for (int ii=0; ii<mm; ii++){
        for (int jj=0; jj<mm; jj++){
        for (int kk=0; kk<mm; kk++){
          real xxc=x0+0.5*ddx+real(ii)*ddx;
          real yyc=y0+0.5*ddy+real(jj)*ddy;
          real zzc=z0+0.5*ddz+real(kk)*ddz;
          real dist=sqrt(pow(xxc,2.0)+pow(yyc,2.0)+pow(zzc-z_init,2.0));
          if (dist>radius){
            itmp=itmp+1;
          }
        }}}
        c[i][j][k]=real(itmp)/real(mm*mm*mm);
      }
    }
    c.bnd_update();
    c.exchange();

    boil::plot->plot(uvw,c,press,"uvw-c-press", 0);
  }

  /*------------------------------+
  |  define static magnetic field |
  +------------------------------*/
  m = Comp::u();
  for_vmijk(B,m,i,j,k){
    B[m][i][j][k] = B0;
  }
  B.exchange_all();

  /* set iint */
  int iint = int(time.current_time()/tint) + 1;
  bool converged_pot=false;

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {

    /*--------------------+
    |  electric potential |
    +--------------------*/
    if (B0!=0 ) {
      real res_pt=1e-4;
      if(!converged_pot){  // reset only when multigrid_pt didn't converged
        boil::oout<<"reset potential\n";
        res_pt=1e-6;
        for_avijk(pot,i,j,k) // in the previous time step
          pot[i][j][k]=0.0;
      }
      pt.discretize();
      pt.coarsen();
      uvw.exchange_all(); // velocity in the corners are used for RHS: u x B
      if(multigrid_pt.vcycle(ResRat(res_pt))) {
        converged_pot=true;
      } else {
        converged_pot=false;
      }
      pot.bnd_update();
      pot.exchange();
      pt.update_j();

      /* shift pot for visualization */
      real pot_ref=0.0;
      if(d.local_i(boil::BW)>=boil::BW &&
         d.local_j(boil::BW)>=boil::BW &&
         d.local_k(boil::BW)>=boil::BW) {
        pot_ref=pot[boil::BW][boil::BW][boil::BW];
      }
      boil::cart.sum_real(&pot_ref);
      for_avijk(pot,i,j,k) pot[i][j][k]=pot[i][j][k]-pot_ref;
    }

    /*--------------------+
    |  momentum equations |
    +--------------------*/
    /* clear body force */
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k] = 0.0;

    /* lorenz force */
    if (B0 != 0.0) {
      pt.force_lorenz(&xyz);
    }

    /* gravity */
    Comp m=Comp::w();
    for_avmijk(xyz,m,i,j,k) {
      xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * mixed.rho(m,i,j,k);
    }

    /* surface tension */
    conc.tension(&xyz, mixed);

    /* solve momentum */
    ns.discretize();
    pr.discretize();
    pr.coarsen();
    ns.new_time_step();
    ns.grad(press);
    ns.solve(ResRat(1e-6));

    /* solve pressure */
    p = 0.0;
    if (multigrid.vcycle(ResRat(1e-3))) OMS(converged);
    p.exchange();
    ns.project(p);
    press += p;
    press.exchange();
    uvw.exchange();
 
    /*-----------------------------+
    | update color function or VOF |
    +-----------------------------*/
    conc.new_time_step();
    conc.advance();
    conc.totalvol();

    /*------------+
    |  floodfill  |
    +------------*/
#ifdef FLOODFILL
    boil::oout<<"main:floodfill\n";
    flood.identify_regions();
#endif

    /* Marino's risign bubble velocity */
    real rising_velocity =0.0;
    real numerator =0.0;
    real denominator =0.0;
    for_vijk(c,i,j,k){
      numerator += (1-c[i][j][k])*(uvw[Comp::w()][i][j][k]+uvw[Comp::w()][i][j][k+1])
                   /2.0*c.dV(i,j,k);
      denominator+=(1-c[i][j][k])*c.dV(i,j,k);
    }
    boil::cart.sum_real(&numerator);
    boil::cart.sum_real(&denominator);
    rising_velocity= numerator/denominator;
    boil::oout<<"rising_velocity= "<<time.current_time()<<" "<<rising_velocity<<"\n";

    /* bubble position using VOF */
    conc.front_minmax();

    if(conc.topo->get_zmaxft()>0.98*LZ){
      boil::oout<<"#main:bubble reached top\n";
      exit(0);
    }

    /*-------------+
    |  dt control  |
    +-------------*/
    real cflmax=ns.cfl_max();
    time.control_dt(cflmax, cfl_limit, dt);

     /*---------------------------+
    |  output for visualization  |
    +---------------------------*/
    if((time.current_time()) / (tint) >= real(iint) ) {
      uvw.exchange_all();
      boil::plot->plot(uvw,c,press,"uvw-c-press", iint);
      boil::plot->plot(xyz,c,"xyz-c",iint);
      if(B0!=0.0){
        boil::plot->plot(J,pot,"J-pot",iint);
        pt.output_j(iint);
      }
      iint++;
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % nint==0) {
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
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      c    .save("c",  time.current_step());
    }

    if( boil::timer.current_min() > (wmin-30)
      || time.current_step()==time.total_steps()) {
      if( boil::cart.iam()==0) {
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
      }
      uvw  .save("uvw",   time.current_step());
      press.save("press", time.current_step());
      c    .save("c",  time.current_step());
      uvw  .rm("uvw",   ts);
      press.rm("press", ts);
      c    .rm("c",  ts);
      boil::timer.stop();
      boil::timer.report();
      exit(0);
    }
  }

  boil::timer.stop();
  boil::timer.report();
}
