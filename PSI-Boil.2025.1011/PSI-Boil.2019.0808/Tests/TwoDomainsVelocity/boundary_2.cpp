  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/

  for_m(m) {
    uvw_2.bc(m).add( BndCnd( Dir::imin(), BndType::insert() ) );
    uvw_2.bc(m).add( BndCnd( Dir::imax(), BndType::outlet() ) );
    uvw_2.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw_2.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_2.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw_2.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  //press_2.bc().add( BndCnd( Dir::imin(), BndType::insert() ) );
  press_2.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  press_2.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  press_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press_2.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press_2.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  p_2.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p_2.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p_2.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  f_2    = p_2.shape();
  q_2    = p_2.shape();
  mdot_2 = p_2.shape();
  mu_t_2 = p_2.shape();

  t_2.bc().add( BndCnd( Dir::imin(), BndType::insert() ) );
  //t_2.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  t_2.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  t_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  t_2.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 110.0 ) );
  t_2.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  c_2.bc().add( BndCnd( Dir::imin(), BndType::insert() ) );
  c_2.bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
  c_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c_2.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c_2.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  step_2 = c_2.shape();

  wd_2.bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
  wd_2.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  wd_2.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  wd_2.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  wd_2.bc().add( BndCnd( Dir::kmin(), BndType::neumann()  ) );
  wd_2.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );

