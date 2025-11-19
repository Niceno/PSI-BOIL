#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  #include "domain-2.cpp"
  const real mu =   0.0003333333333333333333333333333;
	
  Times time(20000, 0.015); /* ndt, dt */
	
  /*--------------------------------+
  |  choose the output file format  |
  +--------------------------------*/
  boil::plot = new PlotGMV();

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d);   // velocity
  Scalar ut(d);  Scalar vt(d);  Scalar wt(d);     // temporary u, v, w
  Scalar u (d);  Scalar v (d);  Scalar w (d);     // averaged  u, v, w
  Scalar uu(d);  Scalar vv(d);  Scalar ww(d);     // averaged  uu,vv,ww
  Scalar uv(d);  Scalar uw(d);  Scalar vw(d);     // averaged  uv,uw,vw

  u  = 0; v  = 0; w  = 0;
  uu = 0; vv = 0; ww = 0;
  uv = 0; uw = 0; vw = 0;

  u.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  u.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );

  /*------------+
  |  time loop  |
  +------------*/
  real count = 0;

  time.first_step(20001);
  for(time.start(); time.end(); time.increase()) {

    /*-------+
    |  load  |
    +-------*/
    if( time.current_step() % 100 == 0 ) {
      boil::oout << "##########################" << boil::endl;
      boil::oout << "# GETTING TIME STEP: " << time.current_step() << boil::endl;
      boil::oout << "#-------------------------" << boil::endl;

      count ++;

      uvw.load("uvw", time.current_step());

      for_vijk(u,i,j,k) {
        ut[i][j][k] = 0.5 * (uvw[0][i][j][k] + uvw[0][i+1][j]  [k]  );
        vt[i][j][k] = 0.5 * (uvw[1][i][j][k] + uvw[1][i]  [j+1][k]  );
        wt[i][j][k] = 0.5 * (uvw[2][i][j][k] + uvw[2][i]  [j]  [k+1]);
      }

      u  += ut;     v  += vt;     w  += wt;
      uu += ut*ut;  vv += vt*vt;  ww += wt*wt;
      uv += ut*vt;  uw += ut*wt;  vw += vt*wt;
    }
  }

  u  /= count;  v  /= count;  w  /= count;
  uu /= count;  vv /= count;  ww /= count;
  uv /= count;  uw /= count;  vw /= count;

  uu -= u*u;  vv -= v*v;  ww -= w*w;
  uv -= u*v;  uv -= u*v;  vw -= v*w;

  boil::oout << "Finished time averaging" << boil::endl;

  /*----------------+
  |  space average  |
  +----------------*/

  real u_ =0.0, v_ =0.0, w_ =0.0; 
  real uu_=0.0, vv_=0.0, ww_=0.0; 
  real uv_=0.0, vw_=0.0, uw_=0.0; 

  for_vj(u,j) {
    u_  = u. average_J(j);
    v_  = v. average_J(j);
    w_  = w. average_J(j);
    uu_ = uu.average_J(j);
    vv_ = vv.average_J(j);
    ww_ = ww.average_J(j);
    uv_ = uv.average_J(j);
    vw_ = vw.average_J(j);
    uw_ = uw.average_J(j);
    for_vik(u,i,k) {
      u [i][j][k] = u_ ;
      v [i][j][k] = v_ ;
      w [i][j][k] = w_ ;
      uu[i][j][k] = uu_;
      vv[i][j][k] = vv_;
      ww[i][j][k] = ww_;
      uv[i][j][k] = uv_;
      vw[i][j][k] = vw_;
      uw[i][j][k] = uw_;
    }
  }

  /*------------+
  |  normalize  |
  +------------*/
  int j; 
  real du_dy, du_dy_1=0, du_dy_2=0;

  if( u.bc().type_here(Dir::jmin(), BndType::wall()) ) 
    du_dy_1 = u.average_J(u.sj()) / u.dys( u.sj() );
  if( u.bc().type_here(Dir::jmax(), BndType::wall()) ) 
    du_dy_2 = u.average_J(u.ej()) / u.dyn( u.ej() );

  boil::cart.sum_real( &du_dy_1 );
  boil::cart.sum_real( &du_dy_2 );

  du_dy = 0.5 * (du_dy_1 + du_dy_2);

  real tau_wall = du_dy * mu;
  real u_tau    = sqrt(tau_wall);

  u /= u_tau;
  uu /= (u_tau * u_tau);
  vv /= (u_tau * u_tau);
  ww /= (u_tau * u_tau);

  /*--------------------+
  |  print out results  |
  +--------------------*/
  Rack r ("r", d, 4, Range<int>(1,128), 68);
  r.print(u);
  r.print(uu);
  r.print(vv);
  r.print(ww);

}
