  /*-----------------+
  |  linear solvers  |
  +-----------------*/
  Krylov * solver_2 = new CG(dom_2);

  /*-----------------+
  |  define solvers  |
  +-----------------*/
  Pressure pr_2( p_2,   f_2,   uvw_2, time, solver_2, &mixed_2 );
  Momentum ns_2( uvw_2, xyz_2,        time, solver_2, &mixed_2 );
  EnthalpyFD en_2(t_2, q_2, c_2, uvw_2, time, solver_2, & mixed_2
                   , tsat, & sapphire_2);
  en_2.convection_set(TimeScheme::forward_euler());
  en_2.diffusion_set(TimeScheme::backward_euler());
  CIPCSL2 conc_2 (c_2,  g_2, uvw_2, time, solver_2);
  conc_2.set_itsharpen(10);
  conc_2.set_globalSharpen();
  Distance di_2(wd_2,  ws_2,  uvw_2, time, solver_2);
  PhaseChange pc_2(mdot_2, t_2, q_2, c_2, g_2, f_2, step_2, uvw_2
                , time, &mixed_2, latent, tsat, &sapphire_2);

  AC sol_2( &pr_2 );
  sol_2.stop_if_diverging(true);
  sol_2.min_cycles(3);
  sol_2.max_cycles(10);

  /*----------------+
  |  wall distance  |
  +----------------*/
  di_2.compute();
  boil::plot->plot(wd_2,"wd_2" );

  /*--------+
  |  model  |
  +--------*/
  Model tm_2;

