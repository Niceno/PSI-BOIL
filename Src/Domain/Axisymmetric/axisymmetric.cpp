//#include "axisymmetric.h"

/******************************************************************************/
Axisymmetric::Axisymmetric(const Grid1D & ogx, const Grid1D & ogz, const real dy,
                           const std::string n,
                           const Decompose dec,
                           const bool print_statistics) :
  Domain(ogx,Grid1D(dy),ogz,n,dec,false),
  ydummy(dy) { 
    grid_y_original = &ydummy;
    check_radial_grid(ogx);

    //boil::oout<<grid_y_local->lx()<<boil::endl;

    if(print_statistics)
      statistics();
}

/******************************************************************************/
Axisymmetric::Axisymmetric(const Grid1D & ogx, const Grid1D & ogz, const real dy,
                           Body * b,
                           const std::string n, 
                           const Decompose dec,
                           const bool print_statistics) :
  Domain(ogx,Grid1D(dy),ogz,b,n,dec,false),
  ydummy(dy) { 
    grid_y_original = &ydummy;
    check_radial_grid(ogx);

    if(print_statistics)
      statistics(body);
}
