  /*------------+
  |  time loop  |
  +------------*/
  bool inertial = time.current_time()<boil::atto;
  
  for(time.start(); time.end(); time.increase()) {

    /* restrict temperature to coarse */
    tpr.restrict_volume_XZ();

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

    /*---------------+
    |  phase change  |
    +---------------*/
    pc_fine.update();

    /* inertial cap */
    if(inertial) {
      real massflux_heat = pc_fine.get_smdot();
      massflux_heat /= conc_fine.topo->get_totarea();
      real massflux_inert = rhov*sqrt(boil::pi/7.*rhov*latent*deltat_nucl
                                      /rhol/tsat0_K);
      boil::oout<<"mflux= "<<massflux_heat<<" "<<massflux_inert<<" "
                           <<massflux_inert/massflux_heat<<boil::endl;
      if(massflux_inert<massflux_heat) {
        mflx.fine *= massflux_inert/massflux_heat;
        mdot.fine *= massflux_inert/massflux_heat;
        pc_fine.sources();
      } else {
        inertial = false;
      }
    }

    /* restrict sources to coarse */
    f.restrict_sum_XZ();
    g.restrict_sum_XZ();

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

    /* we need to get velocity and geometry in the fine space */
    uvw.divergence_free_interpolate_XZ(p,mdot.fine,mixed.fine);

    /*---------------------------+
    |  solve transport equation  |
    +---------------------------*/
    conc_fine.new_time_step();
    conc_coarse.new_time_step();
    conc_coarse.advance_with_extrapolation(false,ResRat(1e-6),uvw.coarse,f.coarse,
                                           &liquid.coarse,&uvw_1);

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
    /* we need to get geometry in the fine space */
    boil::prolongate_color_XZ(conc_coarse,conc_fine);
    conc_fine.color_to_vf();
    conc_fine.ancillary();

    /* in the fine space */
    enthFD_fine.discretize();
    enthFD_fine.new_time_step();
    enthFD_fine.solve(ResRat(1e-16),"enthFD");

    /*-------------+
    |  dt control  |
    +-------------*/
    /* minimum color function */
    conc_coarse.color_minmax();

    /* front */
    conc_fine.front_minmax(Range<real>(0.  ,LX0),
                           Range<real>(-LX0,LX0),
                           Range<real>(0.  ,LZ1));

    time.control_dt(ns.cfl_max(),cfl_limit,dt);

    /*---------------------+
    |  stopping criterion  |
    +---------------------*/
    if(   conc_fine.topo->get_xmaxft()>LX0-dxmin
       || conc_fine.topo->get_zmaxft()>LZ0-dxmin) {
      boil::save_backup(time.current_step(), 1, time,
                        load_scalars, load_scalar_names,
                        load_vectors, load_vector_names);
      boil::rm_backup(ts,
                      load_scalars, load_scalar_names,
                      load_vectors, load_vector_names);

      iint++;
      boil::plot->plot(uvw.coarse,c.coarse,press,
                       "uvwc-cc-press",
                       iint);

      boil::plot->plot(uvw.fine,conc_fine.color(),tpr.fine,mdot.fine,mflx.fine,
                       "uvwf-clrf-tprf-mdot-mflx",
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
      boil::plot->plot(uvw.coarse,c.coarse,press,
                       "uvwc-cc-press",
                       iint);

      boil::plot->plot(uvw.fine,conc_fine.color(),tpr.fine,mdot.fine,mflx.fine,
                       "uvwf-clrf-tprf-mdot-mflx",
                       iint);

      std::fstream output;
      std::stringstream ssp;
      ssp <<"profile-"<<iint<<".txt";
      output.open(ssp.str(), std::ios::out);
      boil::output_profile_xz(conc_coarse.color(),output,Range<int>(NZsol/2+1,NZtot/2));
      boil::cart.barrier();
      output.close();

      std::stringstream ssb;
      ssb <<"bndtpr-"<<iint<<".txt";
      output.open(ssb.str(), std::ios::out);
      boil::output_wall_heat_transfer_xz(tpr.fine,pc_fine.node_tmp(),
                                         solid.fine,output,NXtot);
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
