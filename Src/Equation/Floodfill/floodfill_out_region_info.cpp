#include "floodfill.h"
#include "../../Plot/plot.h"
//#define VERBOSE

/***************************************************************************//**
*  outputs region properties to a file    
*  - center of mass of each region x,y,z
*  - volume of each region
*  - average velocity of each region u,v,w
*******************************************************************************/
void Floodfill::out_region_info() {

  /*----------------------------------------------------+
  |  Output the volume, center of mass, and avg uvw     |
  |  of non hidden regions to a designated output file  |
  +----------------------------------------------------*/
  real odt = time->current_time()-oldtime; //dt between outputs

  if(!boil::cart.iam()) {
#ifdef VERBOSE
    outrgn<<"Info0:t "<<time->current_time()<<" time step "
          <<time->current_step()<<" output period "<<odt
          <<" total number of regions "<<m_vectrgns.size()<<std::endl;
#else
    outrgn<<"Info0:t "<<time->current_time()
          <<" total_regions "<<m_vectrgns.size()<<std::endl;
#endif
  }
  for (int i=0; i<m_vectrgns.size(); i++) {
    Region & r = m_vectrgns[i];
    /* comuvw = dx(com)/dt over dt*50tsteps */
    real cu = (r.x()-r.xold()) / odt;
    real cv = (r.y()-r.yold()) / odt;
    real cw = (r.z()-r.zold()) / odt;
    r.comuvw(cu,cv,cw);
    if (!(r.hiding()) ) { //do not output hidden regions
#ifdef VERBOSE
      outrgn<<"Info1:t "<<time->current_time()<<" ID "<<r.id()
            <<" cvol "<<r.cellvol()
            <<" xyz "<<r.x()<<" "<<r.y()<<" "<<r.z()
            <<" uvw "<<r.u()<<" "<<r.v()<<" "<<r.w()
            <<" comuvw "<<r.comu()<<" "<<r.comv()<<" "<<r.comw()
            <<std::endl;
#else
      outrgn<<"Info1:t "<<time->current_time()<<" ID "<<r.id()
            <<" cvol "<<r.cellvol()
            <<" xyz "<<r.x()<<" "<<r.y()<<" "<<r.z()
            <<std::endl;
#endif
//      <<" hiding "<<r.hiding()<<" tsteps hidden "<<r.get_tsteps_hidden()<<std::endl;
    }
    r.set_oldxyz(); //old xyz com set to current xyz com
  }
  oldtime = time->current_time(); //used to calculate dt between outputs
}
