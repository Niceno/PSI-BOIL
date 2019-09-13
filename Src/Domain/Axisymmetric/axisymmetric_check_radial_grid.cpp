#include "axisymmetric.h"

/***************************************************************************//**
* The grid in the radial direction has to fulfill certain criteria
*******************************************************************************/
void Axisymmetric::check_radial_grid(const Grid1D & gx) {

  real beg = gx.x_min();

  if       (beg< 0.0) {
    boil::oout<<"Radial coordinate can attain only positive values. Fix "
              <<"the radial grid. Exiting."<<boil::endl;
    exit(0);
  } else if(beg==0.0) {
    if(gx.bndgrid1() != BndGrid::wall()) {
       boil::oout<<"Radial grid with zero origin can only have wall "
                 <<"extrapolation at the origin. Exiting."<<boil::endl;
       exit(0);
    }
  } else {
    if(gx.bndgrid1() == BndGrid::symmetry()) {
      boil::oout<<"Radial grid  with non-zero origin cannot be symmetric "
                <<"at the origin. Exiting."<<boil::endl;
      exit(0);
    }
  }

  if(gx.periodic1() == Periodic::yes()||gx.periodicN() == Periodic::yes()) {
       boil::oout<<"Radial grid cannot be periodic. Exiting."<<boil::endl;
       exit(0);
  }

  if(gx.bndgridN() == BndGrid::symmetry()) {
    boil::oout<<"Radial grid cannot be symmetric in the positive radial "
              <<"direction. Exiting."<<boil::endl;
    exit(0);
  }

  if(  gx.bndgrid1() == BndGrid::extrapolate()
     ||gx.bndgridN() == BndGrid::extrapolate()) {
    boil::oout<<"Warning! Constant extrapolation is not well suited for "
              <<"a radial grid!"<<boil::endl;
  }

  boil::oout<<"Axisymmetric domain created, checks for radial grid passed."
            <<boil::endl;

  return;
}
