  /*-------------------+
  |  time-integration  |
  +-------------------*/
  real dt = surftens_dt_coef*pow((vapor.rho()->value()+liquid.rho()->value())
                                 *pow(dxmin,3.0)
                                 /mixed.sigma()->value(),0.5);
  if(case_flag==0) {
    dt = 10.*dxmin;
  }
  boil::oout<<"dxmin= "<<dxmin<<" "<<boil::cart.iam()<<" "<<dt<<"\n";
  Times time(ndt, dt);
  time.set_coef_dec(0.75);
  time.set_dt(dt*initdtcoef);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Krylov * solver = new CG(d, Prec::ic2());
  Krylov * solver_enth = solver;
  //solver_enth = new BiCGS(d, Prec::di());

  /*-------------------+
  |  define equations  |
  +-------------------*/
  /* momentum equation */
  Momentum ns(uvw, xyz, time, solver, &mixed);
  if(mSimple>1)
    ns.convection_set(TimeScheme::backward_euler()); //ns.convection is mandatory
  else
    ns.convection_set(TimeScheme::forward_euler());
  ns.diffusion_set(TimeScheme::backward_euler());

  /* pressure solver */
  Pressure pr(p, f, uvw, time, solver, &mixed);
  AC multigrid( &pr );
  multigrid.stop_if_diverging(multigrid_stop_if_diverging);
  multigrid.use_linf_error(multigrid_use_linf);
  multigrid.min_cycles(multigrid_min_cycles);
  multigrid.max_cycles(multigrid_max_cycles);

  /* color function */
  VOFaxisym conc(c, g, kappa, uvw_1, time, solver);

  conc.set_curv_method(curv_method);
  conc.set_topo_method(topo_method);
  conc.set_use_interp(use_fs_interp);
  conc.set_pressure_extrapolation_parameters(store_pressure_extrap,niter_pressure_extrap);
  conc.set_advection_method(advect_method);

  if(subgrid_method) {
    conc.set_subgrid_method(SubgridMethod::SLICliquid());
  } else {
    conc.set_subgrid_method(SubgridMethod::None());
  }

  conc.set_cangle(cangle);
  conc.set_wall_curv_method(wall_curv_method,Sign::neg(),cangle);
  
  /* supplementary classes for heat transfer */
  TIF tsat(tsat0);

  /* is there conjugate heat transfer? */
  Matter * solid_ptr = NULL;
  if(NZsol>0) {
    solid_ptr = &solid;
  }

  /* function kernel of heat transfer */
  CommonHeatTransfer cht(tpr,conc.topo,tsat,
                         &mixed,solid_ptr);

  if(use_wall_resistance) {
    cht.set_wall_resistance(resistance_liq);
  } else {
    cht.set_int_resistance_liq(resistance_liq);
  }
  if(NZheat>0||NZsol==0) {
  } else {
    cht.set_dirac_wall_source(qflux);
  }

  /* enthalpy equation */
  EnthalpyFDaxisym enthFD(tpr, q, uvw, uvw_1, uvw_2,
                          time, solver_enth, &mixed, cht,
                          solid_ptr);

  enthFD.convection_set(TimeScheme::forward_euler());
  enthFD.diffusion_set(TimeScheme::backward_euler());
  enthFD.set_flux_accuracy_order(ao_efd_conv);

  /* phase change */
  PhaseChange4 pc(mdot, mflx, q,
                  g, f, uvw, cht,
                  time, &mixed, solid_ptr);
  pc.set_accuracy_order(ao_pc);
  pc.set_discard_points_near_interface(discard_points_near_interface);
  pc.set_unconditional_extrapolation(use_unconditional_extrapolation);
