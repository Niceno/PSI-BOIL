/******************************************************************************/
/* ------------ calculated values */
  const int NX0 = int((NXtot-NX_trans*(1.+1./sqrt(AR)))/((LXmult-1.)/AR+1.));
  const int NZ0 = int((NZtot-NZsol-NZ_trans*(1.+1./sqrt(AR)))
                      /((LZmult-1.)/AR+1.));
  const int NZsol0 = int((NZsol-NZsol_trans*(1.+1./sqrt(AR)))
                         /((LZsolmult-1.)/AR+1.));

  boil::oout<<"N0= "<<NX0<<" "<<NZ0<<" "<<NZsol0<<boil::endl;

  const int NX1 = NXtot-NX0-NX_trans;
  const int NZ1 = NZtot-NZ0-NZsol-NZ_trans;
  const int NZsol1 = NZsol-NZsol0-NZsol_trans;

  boil::oout<<"N1= "<<NX1<<" "<<NZ1<<" "<<NZsol1<<boil::endl;

  const real DX0 = NZsol==0 ? LZsol :
                   LZsol/(real(NZsol0)
                          +sqrt(AR)*real(NZsol_trans)+AR*real(NZsol1));
  const real DZ0 = DX0;
  boil::oout<<"DX= "<<DX0<<boil::endl;

  const real LX0 = real(NX0)*DX0; 
  const real LZ0 = real(NZ0)*DZ0; 
  const real LZsol0 = real(NZsol0)*DZ0; 

  boil::oout<<"domain-sizes0= "<<LX0<<" "<<LZ0<<" "<<LZsol0<<boil::endl;

  const int NZheat = std::max(0,int(std::ceil(LZheat/DZ0)));
  const real LZheat_mod = real(NZheat)*DZ0;
  boil::oout<<"NZheat= "<<NZheat<<" "<<LZheat_mod<<boil::endl;

  const real LX_trans = LX0 + sqrt(AR)*DX0*real(NX_trans);
  const real LZ_trans = LZ0 + sqrt(AR)*DZ0*real(NZ_trans);
  const real LZsol_trans = LZsol0 + sqrt(AR)*DZ0*real(NZsol_trans);

  boil::oout<<"domain-sizes_trans= "<<LX_trans<<" "<<LZ_trans<<" "
            <<LZsol_trans<<boil::endl;

  const real LX1 = LX_trans + AR*DX0*real(NX1);
  const real LZ1 = LZ_trans + AR*DZ0*real(NZ1);

  boil::oout<<"domain-sizes1= "<<LX1<<" "<<LZ1<<" "<<LZsol<<boil::endl;

  const real zmax = zmax_mult*LZ1;

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx0( Range<real>(0.0,LX0), NX0, Periodic::no(),
              BndGrid::symmetry(), BndGrid::wall() );
  Grid1D * gxt, * gx1, * gx;
  if(NX_trans>0) {
    gxt = new Grid1D( Range<real>(LX0,LX_trans)
                    , Range<real>(1.0*DX0,AR*DX0)
                    , NX_trans, Periodic::no() );
    if(NX1>0) {
      gx1 = new Grid1D( Range<real>(LX_trans,LX1), NX1, Periodic::no() );

      gx = new Grid1D(gx0, *gxt, *gx1, Periodic::no(), 
                      BndGrid::symmetry(), BndGrid::wall() );
    } else {
      gx = new Grid1D(gx0, *gxt, Periodic::no(), 
                      BndGrid::symmetry(), BndGrid::wall() );
    }
  } else {
    gx = &gx0;
  }

  Grid1D gzf0( Range<real>(0.0,LZ0), NZ0, Periodic::no() );
  Grid1D * gzft, * gzf1, * gzf;
  if(NZ_trans>0) {
    gzft = new Grid1D( Range<real>(LZ0,LZ_trans)
                    , Range<real>(1.0*DZ0,AR*DZ0)
                    , NZ_trans, Periodic::no() );
    if(NZ1>0) {
      gzf1 = new Grid1D( Range<real>(LZ_trans,LZ1), NZ1, Periodic::no() );

      gzf = new Grid1D(gzf0, *gzft, *gzf1, Periodic::no());
    } else {
      gzf = new Grid1D(gzf0, *gzft, Periodic::no());
    }
  } else {
    gzf = &gzf0;
  }

  Grid1D * gzs0, * gzst, * gzs1, * gzs;
  Grid1D * gz;
  if(NZsol>0) {
    gzs0 = new Grid1D( Range<real>(-LZsol0,0.0), NZsol0, Periodic::no() );
    if(NZsol_trans>0) {
      gzst = new Grid1D( Range<real>(-LZsol_trans,-LZsol0)
                      , Range<real>(AR*DZ0,1.0*DZ0)
                      , NZsol_trans, Periodic::no() );
      if(NZsol1>0) {
        gzs1 = new Grid1D( Range<real>(-LZsol,-LZsol_trans),
                           NZsol1, Periodic::no() );

        gzs = new Grid1D(*gzs1, *gzst, *gzs0, Periodic::no());
      } else {
        gzs = new Grid1D(*gzst, *gzs0, Periodic::no());
      }
    } else {
      gzs = gzs0;
    }
    gz = new Grid1D(*gzs, *gzf, Periodic::no());
  } else {
    gz = gzf;
  }

  /*---------+
  |  domain  |
  +---------*/
  Body floor("floor.stl");
  Body * floor_ptr = NULL;
  if(NZsol>0) {
    floor_ptr = &floor;
  }
  TwoLevelAxisymmetric d(*gx,*gz,DX0,floor_ptr);

  const real dxmin = d.coarse().dxyz_min();
  //boil::plot->plot(d.coarse(),"coarsedomain");
