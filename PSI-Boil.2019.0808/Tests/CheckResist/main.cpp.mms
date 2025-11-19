#include "Include/psi-boil.h"

#define CASE 1
/* 
 * 0: linear
 * 1: quadratic
 * 2: sine
 */

/******************************************************************************/
int main(int argc, char * argv[]) {

  boil::timer.start();

  if(argc<3){
    boil::oout<<"Two arguments required!"<<"\n";
    boil::oout<<"./Boil wmin gLevel"<<"\n";
    exit(0);
  }
  const int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const int level=atoi(argv[2]);
  boil::oout<<"level= "<<level<<"\n";

  const real DZ = 0.2/real(level);
  const int NZ = (5+5+5)*level;
  const int NZsol = 5*level;
  
  //const real x0 = 0.77*0.2;//1.03;
  //const real x0 = 0.77*DZ;//1.03;
  const real x0 = 1.03;
  const real Z = 1.;

  //const real qflux = 3.3/(2.*Z);
  const real qflux = 15.3/(2.*Z);
  //const real qflux = 0.0;
  //const real res_wall = 0.2;//2.3;
  const real res_wall = 1.5;
  //const real res_wall = 0.0;//1.5;
  const real res_liq = 0.8;
  //const real res_liq = 0.08;//0.8

#if CASE == 0
  auto analytic = [&](const real t) {
    return x0 - Z*t;
  };

  auto anadot = [&](const real t) {
    return -Z;
  };

  auto anaddot = [&](const real t) {
    return 0.;
  };

  auto anainverse = [&](const real x) {
    return (x0-x)/Z;
  };

  const real tmax = anainverse(0.);
#elif CASE == 1
  auto analytic = [&](const real t) {
    return x0 - Z*t*t;
  };

  auto anadot = [&](const real t) {
    return -2.*Z*t;
  };

  auto anaddot = [&](const real t) {
    return -2.*Z;
  };

  auto anainverse = [&](const real x) {
    return sqrt((x0-x)/Z);
  };

  const real tmax = anainverse(0.);
#elif CASE == 2
  auto analytic = [&](const real t) {
    return x0 - 0.5*Z*sin(2.*boil::pi*t);
  };

  auto anadot = [&](const real t) {
    return -boil::pi*Z*cos(2.*boil::pi*t);
  };

  auto anaddot = [&](const real t) {
    return boil::pi*2.*boil::pi*Z*sin(2.*boil::pi*t);
  };

  auto anainverse = [&](const real x) {
    return asin((x0-x)/Z);
  };

  const real tmax = 1.0;
#endif

  const real t_per_plot = tmax/20.;

  boil::oout<<"tmax= "<<tmax<<" "<<t_per_plot<<boil::endl;

  const real LZ = NZ*DZ;
  const real LZsol = NZsol*DZ;

  const real rhol = 1.;
  const real lambdal = 1.;
  const real cpl = 2.;

  const real rhov = 0.01;
  const real lambdav = 0.01;
  const real cpv = 0.25;

  const real latent = 1.;

  const real rhosol = 4.;
  const real lambdasol = 7.;
  const real cpsol = 5.;

  const real tst0 = 1.;

  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx(DZ);
  Grid1D gz( Range<real>(-LZsol,LZ-LZsol), NZ, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Body * floor_ptr = &floor;
  Domain d(gx, gx, gz,floor_ptr);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Scalar c(d), g(d), kappa(d);        /* concentration */
  Vector uvw(d), xyz(d);              /* velocity */
  Vector uvw_1(d), uvw_2(d);          /* velocity */
  Scalar tpr(d), q(d);                /* temperature */
  Scalar f(d);                        /* pressure src */
  Scalar mdot(d), mflx(d);            /* phase-change rate */

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );

    uvw_1(m) = uvw(m).shape();
    uvw_2(m) = uvw(m).shape();
  }

  c.bc().add( BndCnd( Dir::imin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );

  f     = c.shape();
  g     = c.shape();
  kappa = c.shape();
  q     = c.shape();
  mdot  = c.shape();
  mflx  = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), tst0 ) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
  vapor  .rho   (rhov);
  vapor  .cp    (cpv);
  vapor  .lambda(lambdav);
  liquid.rho   (rhol);
  liquid.cp    (cpl);
  liquid.lambda(lambdal);

  Matter mixed(liquid, vapor,& c);
  mixed.latent(latent);

  Matter solid(d);

  solid.rho    (rhosol);
  solid.cp     (cpsol);
  solid.lambda (lambdasol);

  /*-------------------+
  |  time-integration  |
  +-------------------*/
  const int ndt = 1e7; /* inconsequential */
  const real dxmin = d.dxyz_min();
  const real dt    = 0.1*dxmin;
  const real cfl_limit = 0.2;
  Times time(ndt, dt);
  time.set_coef_dec(0.75);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new BiCGS(d, Prec::di());
  //Krylov * solver = new CG(d, Prec::ic2());

  /*-------------------+
  |  define equations  |
  +-------------------*/
  /* color function */
  VOF conc(c, g, kappa, uvw_1, time, solver);
  conc.set_subgrid_method(SubgridMethod::SLICliquid());

  /* supplementary classes for heat transfer */
  TIF tsat(tst0);

  Matter * solid_ptr = &solid;

  /* function kernel of heat transfer */
  CommonHeatTransfer cht(tpr,conc.topo,tsat,
                         &mixed,solid_ptr);

  if(res_wall>0.) {
    cht.set_wall_resistance(res_wall);
  }
  if(res_liq>0.) {
    cht.set_int_resistance_liq(res_liq);
  }
  if(qflux>0.) {
    cht.set_dirac_wall_source(0.);
  }

  /* enthalpy equation */
  EnthalpyFD enthFD(tpr, q, uvw, uvw_1, uvw_2,
                    time, solver, &mixed, cht,
                    solid_ptr);

  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());

  /* phase change */
  PhaseChange4 pc(mdot, mflx, q,
                  g, f, uvw, cht,
                  time, &mixed, solid_ptr);
  pc.set_accuracy_order(AccuracyOrder::FourthUpwind());
  //pc.set_accuracy_order(AccuracyOrder::FirstUpwind());

/**************** init */
  /* fluid */
  c = 0.;
  for_vijk(c,i,j,k) {
    if(c.zn(k+1)<=x0) {
      c[i][j][k] = 1.0;
    } else if(c.zn(k)<=x0) {
      c[i][j][k] = (x0-c.zn(k))/c.dzc(k);
    }
  }
  c.bnd_update();
  c.exchange_all();
  conc.init();

  /* analytical solution */
  auto tgamma = [&](const real t) {
    return tst0 - res_liq*latent*rhol * anadot(t);
  };

  auto tgammadot = [&](const real t) {
    return -res_liq*latent*rhol * anaddot(t);
  };

  auto M = [&](const real t) {
    return latent*rhol/lambdal/tgamma(t);
  };

  auto Mdot = [&](const real t) {
    return -M(t)/tgamma(t) * tgammadot(t);
  };

  auto E = [&](const real x, const real t) {
    return exp(M(t)*anadot(t)*(x-analytic(t)));
  };

  auto Eprime = [&](const real x, const real t) {
    return M(t)*anadot(t)*E(x,t);
  };

  auto Edoubleprime = [&](const real x, const real t) {
    return M(t)*M(t)*anadot(t)*anadot(t)*E(x,t);
  };

  auto Edot = [&](const real x, const real t) {
    return E(x,t)*( Mdot(t)*(x-analytic(t))*anadot(t)
                   +   M(t)*(x-analytic(t))*anaddot(t)
                   -   M(t)*anadot(t)*anadot(t) );
  };

  auto D = [&](const real t) {
    return -qflux/lambdasol/tgamma(t)/M(t)/E(0.,t);
  };

  auto Ddot = [&](const real t) {
    return -D(t)/tgamma(t)/M(t)/E(0.,t)
           *(tgammadot(t)*M(t)*E(0.,t)
            +tgamma(t)*Mdot(t)*E(0.,t)
            +tgamma(t)*M(t)*Edot(0.,t));
  };

  auto C = [&](const real t) {
    return 1.-lambdal/lambdasol - M(t)*lambdal*res_wall*anadot(t) - D(t);
  };

  auto Cdot = [&](const real t) {
    return - M(t)*lambdal*res_wall*anaddot(t)
           - Mdot(t)*lambdal*res_wall*anadot(t)
           - Ddot(t);
  };

  auto anatprliq = [&](const real x, const real t) {
    return tgamma(t)*E(x,t);
  };

  auto anatprliqprime = [&](const real x, const real t) {
    return tgamma(t)*Eprime(x,t);
  };

  auto anatprliqdoubleprime = [&](const real x, const real t) {
    return tgamma(t)*Edoubleprime(x,t);
  };

  auto anatprliqdot = [&](const real x, const real t) {
    return tgammadot(t)*E(x,t) + tgamma(t)*Edot(x,t);
  };

  auto anatprsol = [&](const real x, const real t) {
    return tgamma(t)*(lambdal/lambdasol+D(t))*E(x,t)
         + tgamma(t)*E(0.,t)*C(t);
  };

  auto anatprsolprime = [&](const real x, const real t) {
    return (lambdal/lambdasol+D(t))*tgamma(t)*Eprime(x,t);
  };

  auto anatprsoldoubleprime = [&](const real x, const real t) {
    return (lambdal/lambdasol+D(t))*tgamma(t)*Edoubleprime(x,t);
  };

  auto anatprsoldot = [&](const real x, const real t) {
    return tgammadot(t)*( (lambdal/lambdasol+D(t))*E(x,t)+E(0.,t)*C(t) )
           +  tgamma(t)*( Ddot(t)*E(x,t)
                        + (lambdal/lambdasol+D(t))*Edot(x,t)
                        + Cdot(t)*E(0.,t)
                        + C(t)*Edot(0.,t) );
  };

  auto anatpr = [&](const real x, const real t) {
     return (x>0.) ? (x>analytic(t)) ? tst0
                                     : anatprliq(x,t)
                   : anatprsol(x,t);
  };

  auto src_liq = [&](const real x, const real t) {
    return cpl*anatprliqdot(x,t) - lambdal*anatprliqdoubleprime(x,t);
  };

  auto src_sol = [&](const real x, const real t) {
    return cpsol*anatprsoldot(x,t) - lambdasol*anatprsoldoubleprime(x,t);
  };

  auto qfunc = [&](const real t) {
      return -qflux*anadot(t);
  };

  /* temperature */
  for_vijk(tpr,i,j,k) {
    tpr[i][j][k] = anatpr(tpr.zc(k),0.);
  }
  tpr.bnd_update();
  tpr.exchange_all();
  
  int iint(0);
  std::fstream output;
  std::stringstream ssp;
  ssp <<"tpr-"<<iint<<".txt";
  output.open(ssp.str(), std::ios::out);
  for_vk(tpr,k) {
    output<<k<<" "<<tpr.zc(k)<<" "<<c[boil::BW][boil::BW][k]<<" "
          <<tpr[boil::BW][boil::BW][k]<<" "
          <<anatpr(tpr.zc(k),0.)
          <<boil::endl;
   }
  boil::cart.barrier();
  output.close();

  conc.front_minmax();
  real zpos = conc.topo->get_zmaxft();

/**************** run */
  for(time.start(); time.end(); time.increase()) {

    real tcurr = time.current_time();

#if 0
    for_vk(tpr,k) {
      boil::oout<<k<<" "<<tpr.zc(k)<<" "<<c[boil::BW][boil::BW][k]<<" "
            <<q[boil::BW][boil::BW][k]<<" "
            <<tpr[boil::BW][boil::BW][k]<<" "
            <<anatpr(tpr.zc(k),tcurr)<<" "
            <<boil::endl;
    }
    boil::oout<<"--- "<<anatprsol(-LZsol,tcurr)<<boil::endl;
#endif

    /* new time step */
    cht.new_time_step();
    conc.new_time_step();

    /*---------------+
    |  phase change  |
    +---------------*/
    pc.update();
    real massflux_heat = pc.get_smdot();
    massflux_heat /= conc.topo->get_totarea();
    real zcalc = massflux_heat/(rhol*latent);
    real difft = zcalc/lambdal;
    boil::oout<<"mflux= "<<time.current_time()<<" "
                         <<massflux_heat<<" "
                         <<zcalc<<" "
                         <<difft<<" "
                         <<-anatprliqprime(analytic(tcurr),tcurr)<<" "
                         <<-anatprliqprime(zpos,tcurr)<<" "
                         <<boil::endl;

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    Comp m = Comp::w();
    boil::oout<<"zpos= "<<time.current_time()<<" "<<zpos<<" "<<analytic(tcurr)<<boil::endl;
    for_avmijk(uvw,m,i,j,k) {
      if(uvw.zc(m,k)>zpos) {
        uvw[m][i][j][k] = massflux_heat*(1./rhov-1./rhol);
      } else {
        uvw[m][i][j][k] = 0.0;
      }
    }
    for_avmijk(uvw,m,i,j,k) {
      uvw_1[m][i][j][k] = 0.0;
      uvw_2[m][i][j][k] = 0.0;//massflux_heat*(1./rhov-1./rhol);
    }

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
#if 1
    conc.advance();
#else
    conc.advance_with_extrapolation(conc.topo->vfold,
                                    true,ResTol(1e-6),uvw,f,
                                    &liquid,&uvw_1,&vapor,&uvw_2);
#endif

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    /* front */
    conc.front_minmax();
    zpos = conc.topo->get_zmaxft();
  
    /* set source */
    for_vijk(q,i,j,k) {
      if(q.zc(k)<=0.) {
        q[i][j][k] = src_sol(q.zc(k),tcurr)*q.dV(i,j,k);
      //} else if(c.zn(k)<=analytic(tcurr)) {
      } else if(q.zc(k)<=zpos) {
        q[i][j][k] = src_liq(q.zc(k),tcurr)*q.dV(i,j,k);
      } else {
        q[i][j][k] = .0;
      }
    }
    q.bnd_update();
    q.exchange_all();

    /* boundary */
    tpr.bc().modify( BndCnd( Dir::kmin(), BndType::dirichlet(), anatprsol(-LZsol,tcurr) ) );
    if(qflux>0.) {
      cht.set_dirac_wall_source(qfunc(tcurr));
    }

    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

    /*-------------+
    |  dt control  |
    +-------------*/
    conc.totalvol();

    /* minimum color function */
    conc.color_minmax();

    real cfl_curr = massflux_heat/rhol*time.dt()/DZ;
    time.control_dt(cfl_curr,cfl_limit,dt);

    /*--------------+
    |  output data  |
    +--------------*/
    bool otpcond = time.current_time() / t_per_plot >= real(iint);
    if(otpcond) {
      iint++;
      std::fstream output;
      std::stringstream ssp;
      ssp <<"tpr-"<<iint<<".txt";
      output.open(ssp.str(), std::ios::out);
      for_vk(tpr,k) {
        output<<k<<" "<<tpr.zc(k)<<" "<<c[boil::BW][boil::BW][k]<<" "
              <<tpr[boil::BW][boil::BW][k]<<" "
              <<anatpr(tpr.zc(k),tcurr)
              <<boil::endl;
      }
      boil::cart.barrier();
      output.close();
    }

    /* terminating condition */
    if(zpos<0.1*DZ||tcurr>tmax)
      break;
  }
  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();

}


