#include "Include/psi-boil.h"
#define USE_VOF

/* boundary conditions */
const real LX = 1.0;
const int  NX = 200;

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*------------------+
  |  plotting format  |
  +------------------*/
  boil::plot = new PlotTEC();

  /*----------+
  |  grid(s)  |
  +----------*/
  //Grid1D gx( Range<real>(-0.5*LX,0.5*LX), NX, Periodic::no() );
  Grid1D gx( Range<real>(0.0,LX), NX, Periodic::no() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gx, gx);
	
  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter fluid(d);

  /*------------------+
  |  define unknowns  |
  +------------------*/
  Vector uvw(d); // velocity
  Scalar c(d), g(d) ; // level set
  Scalar f(d), kappa(d);

  c.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::dirichlet(),0.0 ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::dirichlet(),0.0 ) );

  f = c.shape();


  for_m(m) {
    uvw.bc(m).add( BndCnd( Dir::imin(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::imax(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmin(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::jmax(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmin(), BndType::inlet() ) );
    uvw.bc(m).add( BndCnd( Dir::kmax(), BndType::inlet() ) );
  }

  Krylov * solver = new CG(d, Prec::di());
  int ndt=6000;
  int nint=200;
  Times time(ndt, 0.0005); // ndt, dt


  real pi=acos(-1.0);
#if 0
  /*------------------------+
  |  divergence-free field  |
  +------------------------*/
  real pi=acos(-1.0);
  for(int i=0; i<uvw.ni(); i++)
    for(int j=0; j<uvw.nj(); j++)
      for(int k=0; k<uvw.nk(); k++) {
        real xx,yy,zz;
        xx = uvw.xc(Comp::u(),i);
        yy = uvw.yc(Comp::u(),j);
        zz = uvw.zc(Comp::u(),k);
        uvw[Comp::u()][i][j][k] = 2.0*sin(pi*xx)*sin(pi*xx)
                                * sin(2.0*pi*yy) 
                                * sin(2.0*pi*zz);
        xx = uvw.xc(Comp::v(),i);
        yy = uvw.yc(Comp::v(),j);
        zz = uvw.zc(Comp::v(),k);
        uvw[Comp::v()][i][j][k] = -sin(2.0*pi*xx)
                                * sin(pi*yy)*sin(pi*yy) 
                                * sin(2.0*pi*zz);
        xx = uvw.xc(Comp::w(),i);
        yy = uvw.yc(Comp::w(),j);
        zz = uvw.zc(Comp::w(),k);
        uvw[Comp::w()][i][j][k] = -sin(2.0*pi*xx)
                                * sin(2.0*pi*yy)
                                * sin(pi*zz)*sin(pi*zz);
      }
  uvw.exchange_all(); //set periodic boundary condition.
#endif

  const real radius = 0.15;
  const real xcent = 0.35;
  const real ycent = 0.35;
  const real zcent = 0.35;
  for_vijk(c,i,j,k) {
    real xx=c.xc(i);
    real yy=c.yc(j);
    real zz=c.zc(k);
    real dd = -sqrt((xx-xcent)*(xx-xcent)
                   +(yy-ycent)*(yy-ycent)
                   +(zz-zcent)*(zz-zcent))+radius;
    real eps = LX/NX*1.5;
    if(dd<-eps){
      c[i][j][k]=0.0;
    } else if (dd<eps){
      c[i][j][k]=0.5+dd/(2.0*eps)+1/(2*pi)*sin(pi*dd/eps);
    } else {
      c[i][j][k]=1.0;
    }
  }
  c.exchange_all();
  /*------------------+
  |  define a solver  |
  +------------------*/
  //LevelSet T(c, f, 1.0, 1.0, uvw, time, solver);
  //PhaseField T(c, f, 1.0, 1.0, uvw, time, solver);
#ifdef USE_VOF
  VOF conc  (c, g, kappa, uvw, time, solver);
#else
  CIPCSL2 conc  (c, g, kappa, uvw, time, solver);
#endif
  conc.totalvol();
  boil::plot->plot(uvw,c,"uvw-c", 0);

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

#if 1
    /*------------------------+
    |  divergence-free field  |
    +------------------------*/
    for(int i=0; i<uvw.ni(); i++)
      for(int j=0; j<uvw.nj(); j++)
        for(int k=0; k<uvw.nk(); k++) {
          real const period=3.0;
          real coeft = cos(pi*time.current_time()/period);
          real xx,yy,zz;
          xx = uvw.xc(Comp::u(),i);
          yy = uvw.yc(Comp::u(),j);
          zz = uvw.zc(Comp::u(),k);
          uvw[Comp::u()][i][j][k] = 2.0*sin(pi*xx)*sin(pi*xx)
                                  * sin(2.0*pi*yy)
                                  * sin(2.0*pi*zz) * coeft;
          xx = uvw.xc(Comp::v(),i);
          yy = uvw.yc(Comp::v(),j);
          zz = uvw.zc(Comp::v(),k);
          uvw[Comp::v()][i][j][k] = -sin(2.0*pi*xx)
                                  * sin(pi*yy)*sin(pi*yy)
                                  * sin(2.0*pi*zz) * coeft;
          xx = uvw.xc(Comp::w(),i);
          yy = uvw.yc(Comp::w(),j);
          zz = uvw.zc(Comp::w(),k);
          uvw[Comp::w()][i][j][k] = -sin(2.0*pi*xx)
                                  * sin(2.0*pi*yy)
                                  * sin(pi*zz)*sin(pi*zz) * coeft;
        }
    uvw.exchange_all();
#endif

    conc.advance();
    conc.totalvol();
	  
    if(boil::plot && time.current_step()%nint == 0) {
      boil::plot->plot(uvw,c,"uvw-c", time.current_step());
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-ls.cpp,v 1.17 2009/07/01 14:18:53 niceno Exp $'/
+-----------------------------------------------------------------------------*/
