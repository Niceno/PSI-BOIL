#include "nucleation.h"
using namespace std;

/***************************************************************************//**
*  add nucleation site to Vector sites
*******************************************************************************/
void Nucleation::add(const Site & s) {

  sites.push_back(s);
  set_range(sites);
  vol_area(sites);

  int ns=size()-1;
  int no_daughter = 0;
  int nsd_start = 99999999;
  int nsd_end = 0;
#if 0
  boil::oout<<"vol_area1:ns="<<ns<<" vol_bubble= "<<sites[ns].vol_bubble()
            <<" area_base= "<<sites[ns].area_base()<<"\n";
#endif
#if 0
    std::cout<<"add:x "<<sites[0].x()<<" "<<sites[0].ic()<<" "
             <<sites[0].is()<<" "<<sites[0].ie()<<"\n";
    std::cout<<"add:y "<<sites[0].y()<<" "<<sites[0].jc()<<" "
             <<sites[0].js()<<" "<<sites[0].je()<<"\n";
    std::cout<<"add:z "<<sites[0].z()<<" "<<sites[0].kc()<<" "
             <<sites[0].ks()<<" "<<sites[0].ke()<<"\n";
#endif

  /* add dummy sites for periodic boundary */
  Site stmp=s;
  real xo=stmp.x();
  real yo=stmp.y();
  real zo=stmp.z();
  real LX = clr->domain()->global_max_x() - clr->domain()->global_min_x();
  real LY = clr->domain()->global_max_y() - clr->domain()->global_min_y();
  real act = stmp.active_tpr();
  real zpl = stmp.zplant();

  int nsd;
  /* i-direction */   // crude code
  if (clr->bc().type( Dir::imin(), BndType::periodic())) {
    dsites.push_back( Site(xo-LX, yo, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);

    dsites.push_back( Site(xo+LX, yo, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);
  }
#if 0
  boil::oout<<"vol_area2:ns="<<ns<<" vol_bubble= "<<sites[ns].vol_bubble()
            <<" area_base= "<<sites[ns].area_base()<<"\n";
#endif

  /* j-direction */   // crude code
  if (clr->bc().type( Dir::jmin(), BndType::periodic())) {
    dsites.push_back( Site(xo, yo-LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);

    dsites.push_back( Site(xo, yo+LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);
  }
#if 0
  boil::oout<<"vol_area3:ns="<<ns<<" vol_bubble= "<<sites[ns].vol_bubble()
            <<" area_base= "<<sites[ns].area_base()<<"\n";
#endif

  /* i & j */  // crude code
  if (clr->bc().type( Dir::imin(), BndType::periodic()) &&
      clr->bc().type( Dir::jmin(), BndType::periodic())) {

    dsites.push_back( Site(xo-LX, yo-LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);

    dsites.push_back( Site(xo-LX, yo+LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);

    dsites.push_back( Site(xo+LX, yo-LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);

    dsites.push_back( Site(xo+LX, yo+LY, zo, act, zpl) );
    set_range(dsites);
    nsd = dsites.size()-1;
    vol_area_dummy(ns,nsd);
    dsites[nsd].set_father(ns);
    no_daughter++;
    nsd_start = min(nsd,nsd_start);
    nsd_end = max(nsd,nsd_end);
  }

  /* distribute father's vol and area to daughter */
  for (int nsd=nsd_start; nsd<=nsd_end; nsd++) {
    dsites[nsd].set_vol_bubble(sites[ns].vol_bubble());
    dsites[nsd].set_area_base(sites[ns].area_base());
    //boil::oout<<"add:ns= "<<ns<<" nsd= "<<nsd<<" vol "<<sites[ns].vol_bubble()
    //    <<" "<<dsites[nsd].vol_bubble()<<" area "<<sites[ns].area_base()
    //    <<" "<<dsites[nsd].area_base()<<"\n";
  }

#if 0
  boil::oout<<"nucleation_add:ns="<<ns<<" vol_bubble= "<<sites[ns].vol_bubble()
            <<" area_base= "<<sites[ns].area_base()
	    <<" no_daughter= "<<no_daughter
	    <<" nsd_start= "<<nsd_start<<" nsd_end= "<<nsd_end
	    <<"\n";
#endif
}
