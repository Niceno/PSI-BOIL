  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector xyz(d.coarse());                    /* force */
  Vector uvw_1(d.coarse());                  /* phasic vel */
  Scalar p(d.coarse()), press(d.coarse());   /* pressure */
  Scalar mu_t(d.coarse());                   /* artificial viscosity */

  /* the following variables have to exist on the fine grid: */
  Scalar mdot(d.fine()), mflx(d.fine());     /* phase-change rate */

  /* the following variables have to exist on both levels: */
  TwoLevelVector uvw(d);                     /* velocity */
  TwoLevelScalar c(d), g(d), kappa(d);       /* concentration */
  TwoLevelScalar f(d);                       /* pressure src */
  TwoLevelScalar csol(d);                    /* heater color */
  TwoLevelScalar tpr(d), q(d);               /* temperature */

  /*-----------------------------+
  |  insert boundary conditions  |
  +-----------------------------*/
  for(auto l : uvw.levels) {
    for_m(m) {
      l->bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
      l->bc(m).add( BndCnd( Dir::imax(), BndType::wall() ) );
      l->bc(m).add( BndCnd( Dir::kmin(), BndType::wall() ) );
      l->bc(m).add( BndCnd( Dir::kmax(), BndType::outlet() ) );
      l->bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
      l->bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
    }
  }

  for_m(m)
    uvw_1(m)=uvw.coarse(m).shape();

  for(auto l : c.levels) {
    l->bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    l->bc().add( BndCnd( Dir::imax(), BndType::outlet() ) );
    l->bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
    l->bc().add( BndCnd( Dir::kmax(), BndType::outlet() ) );
    l->bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    l->bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
  }

  press = c.coarse.shape();
  p     = c.coarse.shape();
  mu_t  = c.coarse.shape();

  mdot = c.fine.shape();
  mflx = c.fine.shape();

  for_coarsefine(l) {
    f[l]     = c[l].shape();
    g[l]     = c[l].shape();
    kappa[l] = c[l].shape();
    csol[l]  = c[l].shape();
    q[l]     = c[l].shape();
  }

  for(auto l : tpr.levels) {
    l->bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    l->bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
#ifdef USE_BOTTOM_DIRICHLET
    l->bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(), twall ) );
#else
    l->bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
#endif
    l->bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(),tout) );
    l->bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    l->bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
  }