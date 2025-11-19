  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw_1.bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    uvw_1.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw_1.bc(m).add( BndCnd( Dir::kmax(), BndType::wall() ) );
  }

  press_1.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  press_1.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  press_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  press_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  press_1.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  press_1.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  p_1.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  p_1.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  p_1.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  f_1    = p_1.shape();
  q_1    = p_1.shape();
  mu_t_1 = p_1.shape();

  t_1.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  t_1.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), 110.0 ) );
  t_1.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );

  c_1.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  c_1.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  c_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c_1.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
  c_1.bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
  step_1 = c_1.shape();
  sflag  = c_1.shape();

  wd_1.bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
  wd_1.bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  wd_1.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  wd_1.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  wd_1.bc().add( BndCnd( Dir::kmin(), BndType::neumann()  ) );
  wd_1.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(), 0.0 ) );
  mu_t_1 = wd_1.shape();

