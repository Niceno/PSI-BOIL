#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  /*------------------+
  |  plotting format  |
  +------------------*/
  boil::plot = new PlotTEC(AsNodes::no(), Buffers::yes());

  /*----------+
  |  grid(s)  |
  +----------*/
  Grid1D gx( Range<real>(-1.0,1.0), 10, Periodic::no() );
  Grid1D gy( Range<real>(-0.4,0.4), 4, Periodic::yes() );

  /*---------+
  |  domain  |
  +---------*/
  Domain d(gx, gy, gx);

  /*-------+
  |  time  |
  +-------*/
  int ndt=2;
  int nint=1;
  Times time(ndt, 1.0); // ndt, dt
	
  /*------------------+
  |  define variable  |
  +------------------*/
  Scalar c(d);

  c.bc().add( BndCnd( Dir::imin(), BndType::dirichlet(),
              "(x*x+y*y+z*z)^0.5" ) );
  c.bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  c.bc().add( BndCnd( Dir::imax(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmin(), BndType::wall() ) );
  c.bc().add( BndCnd( Dir::kmax(), BndType::wall() ) );

  /*---------------+
  |  initial cond  |
  +---------------*/
  c.bnd_update();
  c.exchange_all();
  boil::plot->plot(c,"c", 0);
  for_vk(c,k){
    std::cout<<"k= "<<k<<" c[0][1][k]= "<<c[0][1][k]<<"\n";
  }

  for(time.start(); time.end(); time.increase()) {

    boil::oout << "##################" << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME:      " << time.current_time() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "# TIME STEP: " << time.current_step() << boil::endl;
    boil::oout << "#                 " << boil::endl;
    boil::oout << "##################" << boil::endl;

    // set BC as char
    std::string t = std::to_string(time.current_time());
    std::string eq = ("1+(x*x+y*y+z*z)^0.5+");
    std::string eq_t = eq+t;
    std::cout<<eq_t<<"\n";
    char *eqBC = new char[eq_t.length()+1];
    std::strcpy(eqBC, eq_t.c_str());

    // modify BC
    c.bc().modify( BndCnd( Dir::imin(), BndType::dirichlet(), eqBC));

    // update BC
    c.bnd_update();
    c.exchange_all();

    boil::plot->plot(c,"c", time.current_step());
    for_vk(c,k){
      std::cout<<"k= "<<k<<" c[0][1][k]= "<<c[0][1][k]<<"\n";
    }
  }

  boil::oout << "finished" << boil::endl;

  boil::timer.stop();
  boil::timer.report();
}	
/*-----------------------------------------------------------------------------+
 '$Id: main-formula-timeDependBC.cpp,v 1.1 2018/02/20 10:10:17 sato Exp $'/
+-----------------------------------------------------------------------------*/
