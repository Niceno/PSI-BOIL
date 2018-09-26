#include "Include/psi-boil.h"

const real LX =  12.0;
const real LY =   2.0;
const real LZ =   4.0;

/******************************************************************************/
main(int argc, char * argv[]) {

  Times time(40000, 0.04); /* ndt, dt */
	
  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx(Range<real>( 0.0, LX ), 64, Periodic::yes());
  Grid1D gy(Range<real>( -LY/2.0,  LY/2.0), 32, Periodic::no());
  Grid1D gz(Range<real>( 0.0, LZ ), 64, Periodic::yes());

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gz);

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
  Scalar mut(d); Scalar mu(d);                    // temproary and averaged mu_t

  u  = 0.0; v  = 0.0; w  = 0.0;
  uu = 0.0; vv = 0.0; ww = 0.0;
  uv = 0.0; uw = 0.0; vw = 0.0;
  mu = 0.0;

  u.bc().add( BndCnd( Dir::jmin(), BndType::wall() ) );
  u.bc().add( BndCnd( Dir::jmax(), BndType::wall() ) );

  /*------------+
  |  time loop  |
  +------------*/
  real count = 0;

  const Comp U = Comp::u();
  const Comp V = Comp::v();
  const Comp W = Comp::w();

  time.first_step(20001);

  for(time.start(); time.end(); time.increase()) {
    /*-------+
    |  load  |
    +-------*/
    if( time.current_step() % 50 == 0 ) {
      boil::oout << "##########################" << boil::endl;
      boil::oout << "# GETTING TIME STEP: " << time.current_step() << boil::endl;
      boil::oout << "#-------------------------" << boil::endl;

      count += 1.0;

      uvw.load("uvw", time.current_step());

      for_vijk(u,i,j,k) {
        ut[i][j][k] = 0.5 * (uvw[U][i][j][k] + uvw[U][i+1][j]  [k]  );
        vt[i][j][k] = 0.5 * (uvw[V][i][j][k] + uvw[V][i]  [j+1][k]  );
        wt[i][j][k] = 0.5 * (uvw[W][i][j][k] + uvw[W][i]  [j]  [k+1]);
      }

      u  += ut;     v  += vt;     w  += wt;
      uu += ut*ut;  vv += vt*vt;  ww += wt*wt;
      uv += ut*vt;  uw += ut*wt;  vw += vt*wt;

      mut.load("mu_t", time.current_step());
      mu += mut;
    }
  }

  u  /= count;  v  /= count;  w  /= count;
  uu /= count;  vv /= count;  ww /= count;
  uv /= count;  uw /= count;  vw /= count;

  mu /= count;

  uu -= u*u;  vv -= v*v;  ww -= w*w;
  uv -= u*v;  uw -= u*w;  vw -= v*w;

  /*----------------+
  |  space average  |
  +----------------*/

  real u_ =0.0, v_ =0.0, w_ =0.0; 
  real uu_=0.0, vv_=0.0, ww_=0.0; 
  real uv_=0.0, vw_=0.0, uw_=0.0; 
  real mu_=0.0;

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
 
    mu_ = mu.average_J(j); 
    for_vik(mu,i,k) 
      mu[i][j][k] = mu_;
  }

  /*------------+
  |  normalize  |
  +------------*/
  real u_tau = 0.06;

  u /= u_tau;
  uu /= (u_tau * u_tau);
  vv /= (u_tau * u_tau);
  ww /= (u_tau * u_tau);

  for_vijk(u,i,j,k) {
    uu[i][j][k] = sqrt(uu[i][j][k]);
    vv[i][j][k] = sqrt(vv[i][j][k]);
    ww[i][j][k] = sqrt(ww[i][j][k]);
  }

  /*--------------------+
  |  print out results  |
  +--------------------*/
  Rack r ("r", d, 4, Range<int>(1,32), 4);
  r.print(u);
  r.print(uu);
  r.print(vv);
  r.print(ww);
  r.print(mu);
}	

/*-----------------------------------------------------------------------------+
 '$Id: main-y-stat.cpp,v 1.2 2016/02/24 09:36:21 buetikofer_j Exp $'/
+-----------------------------------------------------------------------------*/
