  /*------------+
  |  time loop  |
  +------------*/
  bool inertial = time.current_time()<boil::atto;
  Range<real> ml_range_z(0e-6,6e-6);
  auto cap_frac = [&](const real x, const real z,
                      const Range<real> & xr, const Range<real> zr,
                      const real val) {
    if(z>zr.last()||!xr.exists()) {
      return val;
    } else {
      return val*xr.fraction(x);
    }
  };
  auto cap_val  = [&](const real x, const real z, 
                      const Range<real> & xr, const Range<real> zr,
                      const real val0, const real val1) {
    if(z>zr.last()||!xr.exists()) {
      return val1;
    } else {
      return val0+(val1-val0)*xr.fraction(x);
    }
  };

  conc_coarse.front_minmax(Range<real>(0.  ,LX1),
                           Range<real>(-LX0,LX0),
                           Range<real>(0.  ,LZ1));
  
  for(time.start(); time.end(); time.increase()) {

#ifdef USE_BOTTOM_DIRICHLET
  #include "update_tpr_bnd.cpp"
#endif

    /* test temperature */
    real tsol_max(-boil::unreal), tsol_min(boil::unreal);
    real tliq_max(-boil::unreal), tliq_min(boil::unreal);
    real tvap_max(-boil::unreal), tvap_min(boil::unreal);
    for_vijk(tpr.coarse,i,j,k) {
      real tval = tpr.coarse[i][j][k];
      if(tpr.coarse.domain()->ibody().on(i,j,k)) {
        if(cht_coarse.topo->above_interface(i,j,k)) {
          if(tval>tliq_max)
            tliq_max = tval;
          if(tval<tliq_min)
            tliq_min = tval;
        } else {
          if(tval>tvap_max)
            tvap_max = tval;
          if(tval<tvap_min)
            tvap_min = tval;
        }
      } else {
        if(tval>tsol_max)
          tsol_max = tval;
        if(tval<tsol_min)
          tsol_min = tval;
      }
    }
    boil::cart.max_real(&tsol_max);
    boil::cart.max_real(&tliq_max);
    boil::cart.max_real(&tvap_max);
    boil::cart.min_real(&tsol_min);
    boil::cart.min_real(&tliq_min);
    boil::cart.min_real(&tvap_min);
    boil::oout<<"tprextrema= "<<time.current_time()<<" "
              <<tsol_min<<" "<<tsol_max<<" "
              <<tliq_min<<" "<<tliq_max<<" "
              <<tvap_min<<" "<<tvap_max<<" "
              <<boil::endl;
    if(tliq_min<-0.1||tvap_min<-0.1||tsol_min<-0.1) {
      boil::oout<<"temperature instability. exiting."<<boil::endl;
      iint++;
      boil::plot->plot(uvw.coarse,c.coarse,tpr.coarse,mdot.coarse,mflx.coarse,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);
      exit(0);
    }

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

    /* transition zone */
    real xmaxbub = conc_coarse.topo->get_xmaxft();
    real xmaxml = output_row(ml_range_z.last(),c.coarse,true);
    real xtrml = 1.5*xmaxml - 0.5*xmaxbub;
    boil::oout<<"MLranges= "<<time.current_time()<<" "<<xtrml
                            <<" "<<xmaxml<<" "<<xmaxbub<<boil::endl;
    Range<real> ml_range_x(xtrml,xmaxml);

    real massflow_cap(0.0), massflow_ml(0.0);
    real are_cap(0.0), are_ml(0.0);
    for_vijk(c.coarse,i,j,k) {
      if(conc_coarse.topo->interface(i,j,k)) {
        real a = conc_coarse.topo->get_area(i,j,k);
        real a_cap = cap_frac(c.coarse.xc(i),
                              c.coarse.zc(k),
                              ml_range_x, ml_range_z,
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
    
    real massflux_cap(0.0), massflux_ml(0.0);
    if(are_cap>0.0)
      massflux_cap = massflow_cap/are_cap;
    if(are_ml>0.0)
      massflux_ml = massflow_ml/are_ml;

#if 1
    for_vijk(c.coarse,i,j,k) {
      if(conc_coarse.topo->interface(i,j,k)) {
        mflx.coarse[i][j][k] = cap_val(c.coarse.xc(i),
                                       c.coarse.zc(k),
                                       ml_range_x, ml_range_z,
                                       mflx.coarse[i][j][k],
                                       massflux_cap);
      }
    }
#else
    mflx.coarse = massflux_heat;
#endif
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
    ns.new_time_step(&f.coarse);

    ns.grad(press);
    ns.solve(ResRat(1e-14));

    p = 0.0;
    if(multigrid.cycle(multigrid_cycle0,
                       multigrid_cycle1,
                       multigrid_rt,
                       multigrid_rr,
                       multigrid_mi))
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
    conc_coarse.advance_with_extrapolation(false,ResRat(1e-6),uvw.coarse,f.coarse,
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
