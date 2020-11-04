  /*------------+
  |  time loop  |
  +------------*/
  bool inertial = time.current_time()<boil::atto;
  real z0 = 4e-6;
  real z1 = 7e-6;
  auto cap_frac = [&](const real x, const real z,
                      const real val) {
    return val*std::max(0.0,
                      std::min(1.0,
                        (z-z0)/(z1-z0)));
  };
  
  for(time.start(); time.end(); time.increase()) {

#ifdef USE_BOTTOM_DIRICHLET
  #include "update_tpr_bnd.cpp"
#endif

    /*-------------------+
    |  reset body force  |
    +-------------------*/
    for_m(m)
      for_avmijk(xyz,m,i,j,k)
        xyz[m][i][j][k]=0.0;

    /* gravity force */
    Comp m = Comp::w();
    for_vmijk(xyz,m,i,j,k) {
      real phil=std::max(0.0,std::min(1.0,c.coarse[i][j][k]));
      real phiv=1.0-phil;
      real deltmp=tpr.coarse[i][j][k]-tsat0;
      real rhomix = phil*boil::rho(liquid.coarse.rho()->value(),
                                   liquid.coarse.beta()->value(),deltmp)
                  + phiv*boil::rho(vapor.coarse.rho()->value(),
                                   vapor.coarse.beta()->value(),deltmp);
      if(xyz.domain()->ibody().on(m,i,j,k))
        xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
    }

    /* surface tension */
    conc_coarse.tension(&xyz, mixed.coarse,conc_coarse.color());
    conc_coarse.output_cangle_2d(Comp::i(),Comp::k(),Sign::neg());

    /* boundary temperature */
    cht_coarse.new_time_step();

    /*---------------+
    |  phase change  |
    +---------------*/
    pc_coarse.update();

    real massflow_heat = pc_coarse.get_smdot();
    real massflux_heat = massflow_heat/conc_coarse.topo->get_totarea();
    real massflux_inert = rhov*sqrt(boil::pi/7.*rhov*latent*deltat_nucl
                                    /rhol/tsat0_K);

    real massflow_cap(0.0), massflow_ml(0.0);
    real are_cap(0.0), are_ml(0.0);
    for_vijk(c.coarse,i,j,k) {
      if(conc_coarse.topo->interface(i,j,k)) {
        real a = conc_coarse.topo->get_area(i,j,k);
        real a_cap = cap_frac(c.coarse.xc(i),
                              c.coarse.zc(k),
                              a);
        real a_ml = a-a_cap;  

        massflow_cap += mflx.coarse[i][j][k]*a_cap;
        massflow_ml += mflx.coarse[i][j][k]*a_ml;
        are_cap += a_cap;
        are_ml += a_ml;
      }
    }
    boil::cart.sum_real(&massflow_cap);
    boil::cart.sum_real(&massflow_ml);
    boil::cart.sum_real(&are_cap);
    boil::cart.sum_real(&are_ml);
    
    real massflux_cap = massflow_cap/are_cap;
    real massflux_ml = massflow_ml/are_ml;

    mflx.coarse = massflux_heat;
    mflx.coarse.bnd_update();
    mflx.coarse.exchange();
    pc_coarse.finalize();

    boil::oout<<"mflux= "<<time.current_time()<<" "
                         <<massflux_heat<<" "<<massflux_inert<<" "
                         <<massflux_inert/massflux_heat<<" | "
                         <<massflow_heat<<" "<<massflow_cap<<" "<<massflow_ml<<" | "
                         <<massflux_cap<<" "<<massflux_ml
                         <<boil::endl;

    ns.vol_phase_change(&f.coarse);

    /*--------------------------+
    |  solve momentum equation  |
    +--------------------------*/
    /* essential for moving front */
    ns.discretize();
    pr.discretize();
    pr.coarsen();

    /* momentum */
    ns.new_time_step();

    ns.grad(press);
    ns.solve(ResRat(1e-14));

    p = 0.0;
    if(multigrid.cycle(multigrid_cycle0,
                       multigrid_cycle1,
                       multigrid_rr,multigrid_mi))
      OMS(converged);

    p.exchange();
    ns.project(p);
    press += p;

    /* shift pressure */
    real pmin=1.0e+300;
    for_vijk(press,i,j,k) {
      if(press.domain()->ibody().on(i,j,k)) {
        if(pmin>press[i][j][k])
          pmin=press[i][j][k];
      }
    }
    boil::cart.min_real(&pmin);

    for_vijk(press,i,j,k) {
      if(press.domain()->ibody().on(i,j,k)) {
        press[i][j][k] -= pmin;
      } else {
        press[i][j][k] = 0.0;
      }
    }

    press.bnd_update();
    press.exchange_all();

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc_coarse.new_time_step();
    conc_coarse.advance_with_extrapolation(false,ResRat(1e-9),uvw.coarse,f.coarse,
                                           &liquid.coarse,&uvw_1,&vapor.coarse,&uvw_2);

    for_avk(c.coarse,k) {
      if(c.coarse.zc(k)>=(zmax-c.coarse.dzc(k))) {
        for_avij(c.coarse,i,j) {
          c.coarse[i][j][k]= 1.0;
        }
      }
    }

    c.coarse.bnd_update();
    c.coarse.exchange_all();
    conc_coarse.ancillary();
    conc_coarse.totalvol();

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    tprold.coarse = tpr.coarse;
#if 1
    enthFD_coarse.discretize();
    enthFD_coarse.new_time_step();
    enthFD_coarse.solve(ResRat(1e-16),"enthFD");
#else
    enthFD_coarse.discretize();
    //enthFD_coarse.new_time_step();
    for(int i(0); i<3;++i) {
      for(int j(0); j<3;++j) {
        enthFD_coarse.convective_time_step(tprold.coarse);
        if(j==0) {
          tprap1.coarse = tpr.coarse;
        } else {
          tpr.coarse = (tprap1.coarse+tpr.coarse);
          tpr.coarse /= 2.;
        }
      }
      enthFD_coarse.inertial(tpr.coarse,false,Old::no);
      enthFD_coarse.solve(ResRat(1e-16));
      if(i==0) {
        tprap2.coarse = tpr.coarse;
      } else {
        tpr.coarse = (tprap2.coarse+tpr.coarse);
        tpr.coarse /= 2.;
      }
    }
#endif

    /*-------------+
    |  dt control  |
    +-------------*/
    /* minimum color function */
    conc_coarse.color_minmax();

    /* front */
    conc_coarse.front_minmax(Range<real>(0.  ,LX1),
                           Range<real>(-LX0,LX0),
                           Range<real>(0.  ,LZ1));

    time.control_dt(ns.cfl_max(),cfl_limit,dt);

    /*---------------------+
    |  stopping criterion  |
    +---------------------*/
    if(   conc_coarse.topo->get_xmaxft()>LX0-dxmin
       || conc_coarse.topo->get_zmaxft()>LZ0-dxmin) {
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names);

      iint++;
      boil::plot->plot(uvw.coarse,c.coarse,tpr.coarse,mdot.coarse,mflx.coarse,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);

      /* cell-center velocities */
      Scalar u(d.coarse()), v(d.coarse()), w(d.coarse());
      boil::cell_center_velocities(uvw.coarse,u,v,w);
      boil::save_backup(time.current_step(), 1, time,
                        {&u,&v,&w}, {"u","v","w"});

      break;
    }

    /*--------------+
    |  output data  |
    +--------------*/
    bool otpcond = time.current_time() / t_per_plot >= real(iint);
    if(otpcond) {
      iint++;
      boil::plot->plot(uvw.coarse,c.coarse,tpr.coarse,mdot.coarse,mflx.coarse,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);

      std::fstream output;
      std::stringstream ssp;
      ssp <<"profile-"<<iint<<".txt";
      output.open(ssp.str(), std::ios::out);
      boil::output_profile_xz(conc_coarse.color(),output,Range<int>(NZsol/2+1,NZtot/2),
                              Range<int>(-1,-2),LX1);
      boil::cart.barrier();
      output.close();

      real ml_thickness_max = 10e-6;
      std::stringstream ssm;
      ssm <<"microlayer-"<<iint<<".txt";
      output.open(ssm.str(), std::ios::out);
      boil::output_profile_zx(conc_coarse.color(),output,
                              Range<int>(1,NXtot/2),
                              Range<int>(NZsol/2+1,
                                         NZsol/2+ceil(ml_thickness_max/2./DX0)
                                        )
                             );
      boil::cart.barrier();
      output.close();

      std::stringstream ssb;
      ssb <<"bndtpr-"<<iint<<".txt";
      output.open(ssb.str(), std::ios::out);
      if(NZsol>0) {
        boil::output_wall_heat_transfer_xz(cht_coarse,output,NXtot/2);
      } else {
        boil::output_wall_heat_transfer_xz(tpr.coarse,*(conc_coarse.topo),pc_coarse,
                                           output,NXtot/2);
      }
      boil::cart.barrier();
      output.close();
    }

    /*--------------+
    |  backup data  |
    +--------------*/
    if(time.current_step() % n_per_backup == 0) {
      boil::save_backup(time.current_step(), 0, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
    }

    if(boil::timer.current_min() > (wmin-30.0)) {
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names);

      boil::set_irun(0);

      break;
    }
  }
