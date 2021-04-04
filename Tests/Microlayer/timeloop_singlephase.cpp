  /*---------------------------+
  |  solve transport equation  |
  +---------------------------*/
  conc.new_time_step();

  /*------------+
  |  time loop  |
  +------------*/
  for(time.start(); time.end(); time.increase()) {
    cht.new_time_step();

#ifdef USE_MOMEMTUM_SINGLE_PHASE
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
                       multigrid_rt,
                       multigrid_rr,
                       multigrid_mi,
                       multigrid_mstale))
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
#endif

    /*------------------------+
    |  solve energy equation  |
    +------------------------*/
    enthFD.discretize();
    enthFD.new_time_step();
    enthFD.solve(ResRat(1e-16),"enthFD");

    /*-------------+
    |  dt control  |
    +-------------*/
    time.control_dt(ns.cfl_max(),cfl_limit,dt);

    /*---------------------+
    |  stopping criterion  |
    +---------------------*/
    real tprtest(0.);
    pc.update();
    for_vmijk(cht.node_tmp_flu(),Comp::k(),i,j,k) {
      if(fabs(cht.node_tmp_flu().zc(Comp::k(),k))<boil::atto) {
        if(cht.node_tmp_flu().xc(Comp::k(),i)<dxmin) {
          tprtest = cht.node_tmp_flu()[Comp::k()][i][j][k];
        }
      }
    }
    boil::cart.max_real(&tprtest);
    boil::oout<<"tprnucl= "<<time.current_time()<<" "<<tprtest<<boil::endl;

    if(tprtest>tnucl) {
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
      boil::plot->plot(uvw,c,tpr,press,
                       "uvw-c-tpr-press",
                       iint);

      /* output temperature field */
      if(!boil::cart.iam())
        output_to_file(tpr,cht.node_tmp_flu());

      std::fstream output;
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

      break;
    }

    /*--------------+
    |  output data  |
    +--------------*/
    bool otpcond = time.current_time() / t_per_plot >= real(iint);
    if(otpcond) {
      iint++;
      boil::plot->plot(uvw,c,tpr,press,
                       "uvw-c-tpr-press",
                       iint);

      std::fstream output;
      std::stringstream ssb;
      ssb <<"bndtpr-"<<iint<<".txt";
      output.open(ssb.str(), std::ios::out);
      pc.update();
      boil::output_wall_heat_transfer_xz(cht,output,NXtot);
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
