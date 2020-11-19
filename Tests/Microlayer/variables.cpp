  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector xyz(d);                      /* force */
  Vector uvw_1(d), uvw_2(d);          /* phasic vel */
  Scalar p(d), press(d);              /* pressure */
  Scalar mu_t(d);                     /* artificial viscosity, unused */

  Vector uvw(d), uvw_old(d);          /* velocity */
  Scalar c(d), g(d), kappa(d);        /* concentration */
  Scalar f(d);                        /* pressure src */
  Scalar csub(d);                     /* heater color */
  Scalar tpr(d), tprold(d), q(d);     /* temperature */
  Scalar mdot(d), mflx(d);            /* phase-change rate */

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

    uvw_1(m)=uvw(m).shape();
    uvw_2(m)=uvw(m).shape();
    uvw_old(m)=uvw(m).shape();
  }

  c.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

  press = c.shape();
  p     = c.shape();
  mu_t  = c.shape();
  f     = c.shape();
  g     = c.shape();
  kappa = c.shape();
  csub  = c.shape();
  q     = c.shape();
  mdot  = c.shape();
  mflx  = c.shape();

  tpr.bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
  tpr.bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
#ifdef USE_BOTTOM_DIRICHLET
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), twall ) );
#else
  tpr.bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
#endif
  //tpr.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(),tout) );
  tpr.bc().add( BndCnd( Dir::kmax(), BndType::outlet()) );
  tpr.bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
  tpr.bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );

  tprold = tpr.shape();
