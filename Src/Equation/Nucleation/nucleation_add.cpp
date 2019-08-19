#include "nucleation.h"

/***************************************************************************//**
*  add nucleation site to Vector sites
*******************************************************************************/
void Nucleation::add(const Site & s) {

  sites.push_back(s);
  set_range(sites);
#if 0
    std::cout<<"add:x "<<sites[0].x()<<" "<<sites[0].ic()<<" "
             <<sites[0].is()<<" "<<sites[0].ie()<<"\n";
    std::cout<<"add:y "<<sites[0].y()<<" "<<sites[0].jc()<<" "
             <<sites[0].js()<<" "<<sites[0].je()<<"\n";
    std::cout<<"add:z "<<sites[0].z()<<" "<<sites[0].kc()<<" "
             <<sites[0].ks()<<" "<<sites[0].ke()<<"\n";
#endif

  /* add dummy sites for periodic boundary */
  Site stmp = s;
  real xo=stmp.x();
  real yo=stmp.y();
  real zo=stmp.z();
  real LX = clr->domain()->global_max_x() - clr->domain()->global_min_x();
  real LY = clr->domain()->global_max_y() - clr->domain()->global_min_y();
  real act = stmp.active_tpr();
  real zpl = stmp.zplant();

  int  sID=size()-1;
  /* i-direction */   // crude code
  if (clr->bc().type( Dir::imin(), BndType::periodic())) {
    dummy_add( Site(xo-LX, yo, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
#if 0
    std::cout<<"add:x"<<dsites[0].x()<<" "<<dsites[0].ic()<<" "
             <<dsites[0].is()<<" "<<dsites[0].ie()<<"\n";
    std::cout<<"add:y"<<dsites[0].y()<<" "<<dsites[0].jc()<<" "
             <<dsites[0].js()<<" "<<dsites[0].je()<<"\n";
    std::cout<<"add:z"<<dsites[0].z()<<" "<<dsites[0].kc()<<" "
             <<dsites[0].ks()<<" "<<dsites[0].ke()<<"\n";
#endif
    dummy_add( Site(xo+LX, yo, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
  }

  /* j-direction */   // crude code
  if (clr->bc().type( Dir::jmin(), BndType::periodic())) {
    dummy_add( Site(xo, yo-LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
    dummy_add( Site(xo, yo+LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
  }

  /* i & j */  // crude code
  if (clr->bc().type( Dir::imin(), BndType::periodic()) &&
      clr->bc().type( Dir::jmin(), BndType::periodic())) {
    dummy_add( Site(xo-LX, yo-LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
    dummy_add( Site(xo-LX, yo+LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
    dummy_add( Site(xo+LX, yo-LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
    dummy_add( Site(xo+LX, yo+LY, zo, act, zpl) );
    set_range(dsites);
    dsites[dsites.size()-1].set_father(sID);
  }

}
