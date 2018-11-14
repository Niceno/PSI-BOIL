  /*-----------------+
  |  linear solvers  |
  +-----------------*/
  Krylov * solver_1 = new CG(dom_1);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr_1( p_1,   f_1,   uvw_1, time, solver_1, &mixed_1 );
  Momentum ns_1( uvw_1, xyz_1,        time, solver_1, &mixed_1 );
  EnthalpyFD en_1(t_1, q_1, c_1, uvw_1, time, solver_1, & mixed_1
                   , tsat, & sapphire_1);
  en_1.convection_set(TimeScheme::forward_euler());
  en_1.diffusion_set(TimeScheme::backward_euler());
  CIPCSL2 conc_1 (c_1,  g_1, uvw_1, time, solver_1);
  conc_1.set_itsharpen(10);
  conc_1.set_globalSharpen();
  Distance di_1(wd_1,  ws_1,  uvw_1, time, solver_1);
  //PhaseChange pc_1(mdot_1, t_1, q_1, c_1, g_1, f_1, step_1, uvw_1
  //              , time, &mixed_1, latent, tsat, &sapphire_1);

  AC sol_1( &pr_1 );
  sol_1.stop_if_diverging(true);
  sol_1.min_cycles(3);
  sol_1.max_cycles(10);

  /*----------------+
  |  wall distance  |
  +----------------*/
  di_1.compute();
  boil::plot->plot(wd_1,"wd_1" );

  /*--------+
  |  model  |
  +--------*/
  Model tm_1;

