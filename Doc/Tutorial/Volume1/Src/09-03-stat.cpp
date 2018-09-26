#include "Include/psi-boil.h"

#include <vector>

/****************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*------------------------------+
  |  grids, obstacles and domain  |
  +------------------------------*/
  #include "09-03-common.h"
	
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

  /* timer */
  Times time(80000, 0.00002); /* ndt, dt */
  time.first_step(20000);

  /*------------+
  |  time loop  |
  +------------*/
  real count = 0;

  const Comp U = Comp::u();
  const Comp V = Comp::v();
  const Comp W = Comp::w();
 
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
        ut[i][j][k] = 0.5 * (uvw[U][i][j][k] + uvw[U][i+1][j]  [k]  );
        vt[i][j][k] = 0.5 * (uvw[V][i][j][k] + uvw[V][i]  [j+1][k]  );
        wt[i][j][k] = 0.5 * (uvw[W][i][j][k] + uvw[W][i]  [j]  [k+1]);
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

  /*-------+
  |  plot  |
  +-------*/
  boil::plot->plot(u, v, w,   "velocity-mean", time.current_step()-1);
  boil::plot->plot(uu,vv,ww,  "stresses-mean", time.current_step()-1);

  boil::timer.stop();
  boil::timer.report();
}

/*-----------------------------------------------------------------------------+
 '$Id: 09-03-stat.cpp,v 1.5 2009/07/13 09:29:32 niceno Exp $'/
+-----------------------------------------------------------------------------*/
