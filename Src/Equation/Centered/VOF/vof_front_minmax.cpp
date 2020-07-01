#include "vof.h"

/******************************************************************************/
void VOF::front_minmax() {
  front_minmax(Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa),
               Range<real>(-boil::exa, boil::exa));

  return;
}

/******************************************************************************/
void VOF::front_minmax(Range<real> xr,
                       Range<real> yr,
                       Range<real> zr ) {
  topo->front_minmax(xr,yr,zr);

  std::cout.setf(std::ios_base::scientific);
  boil::oout<<"x-min-front= "<<time->current_time()<<" "
            <<topo->get_xminft()<<" "<<topo->get_xmaxft()<<" "
            <<topo->get_yminft()<<" "<<topo->get_ymaxft()<<" "
            <<topo->get_zminft()<<" "<<topo->get_zmaxft()<<" "
            <<"\n";
  std::cout.unsetf(std::ios_base::floatfield);

  return;
}
