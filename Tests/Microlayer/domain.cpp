/******************************************************************************/
/* ------------ calculated values */
  const real DX0 = LZsol/real(NZsol);
  boil::oout<<"DX= "<<DX0<<boil::endl;
  const real DZ0 = DX0;

  const int NZheat = std::max(0,int(std::ceil(LZheat/DZ0)));
  const real LZheat_mod = real(NZheat)*DZ0;
  boil::oout<<"NZheat= "<<NZheat<<" "<<LZheat_mod<<boil::endl;

  const int NX0 = finefact*NXtot;
  const int NZ0 = finefact*(NZtot-NZsol);

  const int NX1 = NXtot-NX0;
  const int NZ1 = NZtot-NZ0-NZsol;

  const real LX0 = real(NX0)*DX0; 
  const real LZ0 = real(NZ0)*DZ0; 

  const real LX1 = 1.5*LX0;
  const real LZ1 = 1.5*LZ0;

  const real zmax = zmax_mult*LZ1;

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx0( Range<real>(0.0,LX0), NX0, Periodic::no() );
  Grid1D gx1( Range<real>(LX0,LX1)
            , Range<real>(1.1*DX0,3.1*DX0)
            , NX1, Periodic::no() );

  Grid1D gx(gx0, gx1, Periodic::no(), BndGrid::symmetry(), BndGrid::wall() );

  Grid1D gz0( Range<real>(0.0, LZ0), NZ0, Periodic::no() );
  Grid1D gz1( Range<real>(LZ0,LZ1)
            , Range<real>(1.1*DZ0,3.1*DZ0)
            , NZ1, Periodic::no() );

  Grid1D gzf(gz0, gz1, Periodic::no() );
  Grid1D gzs( Range<real>(-LZsol, 0.0), NZsol, Periodic::no() );

  Grid1D gz(gzs, gzf, Periodic::no());

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  TwoLevelAxisymmetric d(gx,gz,DX0,&floor);

  const real dxmin = d.coarse().dxyz_min();
  boil::plot->plot(d.coarse(),"coarsedomain");