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

  conc.front_minmax(Range<real>(0.  ,LX1),
                    Range<real>(-LX0,LX0),
                    Range<real>(0.  ,LZ1));
  
  for(time.start(); time.end(); time.increase()) {

#ifdef USE_BOTTOM_DIRICHLET
  #include "update_tpr_bnd.cpp"
#endif

#include "tpr_extrema.cpp"

#ifdef USE_DYNCA
  #include "dynamic_ca.cpp"
#endif

/**********************************************************************/

    /* store velocity & temperature */
    if(mSimple>1) {
      for_m(m)
        uvw_old(m) = uvw(m);
      tprold = tpr;
    } 

    /* new time step */
    cht.new_time_step();
    conc.new_time_step();
    conc.output_cangle_2d(Comp::i(),Comp::k(),Sign::neg());

    /* inner loop */
    for(int mloop=0; mloop<mSimple; mloop++) {
      /*-------------------+
      |  reset body force  |
      +-------------------*/
      for_m(m)
        for_avmijk(xyz,m,i,j,k)
          xyz[m][i][j][k]=0.0;

      /* gravity force */
      Comp m = Comp::w();
      for_vmijk(xyz,m,i,j,k) {
        real phil=std::max(0.0,std::min(1.0,c[i][j][k]));
        real phiv=1.0-phil;
        real deltmp=tpr[i][j][k]-tsat0;
        real rhomix = phil*boil::rho(liquid.rho()->value(),
                                     liquid.beta()->value(),deltmp)
                    + phiv*boil::rho(vapor.rho()->value(),
                                     vapor.beta()->value(),deltmp);
        if(xyz.domain()->ibody().on(m,i,j,k))
          xyz[m][i][j][k] += -gravity * xyz.dV(m,i,j,k) * rhomix;
      }

      /* surface tension */
      conc.tension(&xyz, mixed,conc.color());

      /*---------------+
      |  phase change  |
      +---------------*/
      pc.update();

      /* initial averaging */
      if(conc.topo->get_xmaxft()<150e-6) {
        real massflow_heat = pc.get_smdot();
        real massflux_heat = massflow_heat/conc.topo->get_totarea();

        mflx = massflux_heat;
        mflx.bnd_update();
        mflx.exchange();
        pc.finalize();
      }

      ns.vol_phase_change(&f);

      /*--------------------------+
      |  solve momentum equation  |
      +--------------------------*/
      ns.discretize();
      pr.discretize();
      pr.coarsen();
  
      /* momentum */
      if(mSimple>1) {
        ns.new_time_step(uvw_old,&f);
        ns.convection();
      } else {
        ns.new_time_step(&f);
      }

      ns.grad(press);
      ns.solve(ResRat(1e-14));

      /* pressure */
      p = 0.0;
      if(multigrid.cycle(multigrid_cycle0,
                         multigrid_cycle1,
                         multigrid_rt,
                         multigrid_rr,
                         multigrid_mi,
                         multigrid_mstale))
        OMS(converged);

      p.exchange();
      ns.project(p);
      press += p;

      /*---------------------------+
      |  solve transport equation  |
      +---------------------------*/
      conc.advance_with_extrapolation(conc.topo->vfold,
                                      false,ResTol(1e-6),uvw,f,
                                      &liquid,&uvw_1,&vapor,&uvw_2);

      for_avk(c,k) {
        if(c.zc(k)>=(zmax-c.dzc(k))) {
          for_avij(c,i,j) {
            c[i][j][k]= 1.0;
          }
        }
      }
      c.bnd_update();
      c.exchange_all();
      conc.ancillary();

      /*------------------------+
      |  solve energy equation  |
      +------------------------*/
      if(mloop>0)
        tpr = tprold;
      enthFD.discretize();
      enthFD.new_time_step();
      enthFD.solve(ResRat(1e-16),"enthFD");
    }

#include "shift_pressure.cpp"

#include "eval_massflux.cpp"

    /*-------------+
    |  dt control  |
    +-------------*/
    conc.totalvol();

    /* minimum color function */
    conc.color_minmax();

    /* front */
    conc.front_minmax(Range<real>(0.  ,LX1),
                      Range<real>(-LX0,LX0),
                      Range<real>(0.  ,LZ1));

    time.control_dt(ns.cfl_max(),cfl_limit,dt);
#if 0
    real cap_ts = surftens_dt_coef*conc.topo->capillary_ts(mixed,uvw_1,uvw_2);
    boil::oout<<"cap_ts= "<<time.current_time()<<" "
              <<cap_ts<<" "<<time.dt()<<boil::endl;
    time.control_dt(time.dt(),cap_ts,dt);
#endif

    /*---------------------+
    |  stopping criterion  |
    +---------------------*/
    if(   conc.topo->get_xmaxft()>LX0-dxmin
       || conc.topo->get_zmaxft()>LZ0-dxmin) {
      /* cell-center velocities */
      Scalar u(d), v(d), w(d);
      boil::cell_center_velocities(uvw,u,v,w);
      boil::save_backup(time.current_step(), 1, time,
                        {&u,&v,&w}, {"u","v","w"});

      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names,
                        load_nucls,   load_nucl_names,
                        load_cipcsls, load_cipcsl_names,
                        store_reals,  store_ints);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names,
                      load_nucls,   load_nucl_names,
                      load_cipcsls, load_cipcsl_names,
                      store_reals,  store_ints);

      iint++;
      boil::plot->plot(uvw,c,tpr,mdot,mflx,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);

      break;
    }

    /*--------------+
    |  output data  |
    +--------------*/
    bool otpcond = time.current_time() / t_per_plot >= real(iint);
    if(otpcond) {
      iint++;
      boil::plot->plot(uvw,c,tpr,mdot,mflx,press,
                       "uvw-c-tpr-mdot-mflx-press",
                       iint);

      std::fstream output;
      std::stringstream ssp;
      ssp <<"profile-"<<iint<<".txt";
      output.open(ssp.str(), std::ios::out);
      boil::output_profile_xz(conc.color(),output,Range<int>(NZsol+1,NZtot),
                              Range<int>(-1,-2),LX1);
      boil::cart.barrier();
      output.close();

      real ml_thickness_max = 10e-6;
      std::stringstream ssm;
      ssm <<"microlayer-"<<iint<<".txt";
      output.open(ssm.str(), std::ios::out);
      boil::output_profile_zx(conc.color(),output,
                              Range<int>(1,NXtot),
                              Range<int>(NZsol+1,
                                         NZsol+ceil(ml_thickness_max/DX0)
                                        )
                             );
      boil::cart.barrier();
      output.close();

      std::stringstream ssb;
      ssb <<"bndtpr-"<<iint<<".txt";
      output.open(ssb.str(), std::ios::out);
      if(NZsol>0) {
        boil::output_wall_heat_transfer_xz(cht,output,NXtot);
      } else {
        boil::output_wall_heat_transfer_xz(tpr,*(conc.topo),pc,
                                           output,NXtot);
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
                        load_vectors, load_vector_names,
                        load_nucls,   load_nucl_names,
                        load_cipcsls, load_cipcsl_names,
                        store_reals,  store_ints);
    }

    if(boil::timer.current_min() > (wmin-30.0)) {
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names,
                        load_nucls,   load_nucl_names,
                        load_cipcsls, load_cipcsl_names,
                        store_reals,  store_ints);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names,
                      load_nucls,   load_nucl_names,
                      load_cipcsls, load_cipcsl_names,
                      store_reals,  store_ints);

      boil::set_irun(0);

      break;
    }
  }
