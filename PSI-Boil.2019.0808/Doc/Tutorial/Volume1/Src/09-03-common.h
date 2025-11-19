  /* parameters */
  const real LX  =  0.05999;
  const real LY  =  0.05099;
  const real rho =  1.205;
  const real nu  = 15.11e-6;
  const real mu  = rho * nu; 

  const int NX = 128;
  const int NY = 108;

  const real b_des = 3.86;

  boil::plot = new PlotTEC();

  /*-------------------------------------+
  |  1. grids, immersed body and domain  |
  +-------------------------------------*/
  Grid1D gx(Range<real>( -LX/2.0,  LX/2.0 ), NX, Periodic::yes());
  Grid1D gy(Range<real>(  0,       LY     ), NY, Periodic::no());

  Body cube("09-03-cube.stl");

  Domain d(gx, gx, gy, & cube);
