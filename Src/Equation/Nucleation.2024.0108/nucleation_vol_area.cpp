#include "nucleation.h"
//#define DEBUG

/***************************************************************************//**
* calculate volume and bottom_area of seed bubble
*******************************************************************************/
void Nucleation::vol_area(std::vector<Site> & s) {

  int ns=s.size()-1;

  real vol=0.0;   // volume of seed
  real area_base=0.0;  // area of bubble-base
  
  if (s[ns].contain_range() ) {

#if 0
    if (ns==349) {
    std::cout<<"vol_area:proc="<<boil::cart.iam()
             <<" isie= "<<s[ns].is()<<" "<<s[ns].ie()
             <<" jsje= "<<s[ns].js()<<" "<<s[ns].je()
             <<" kske= "<<s[ns].ks()<<" "<<s[ns].ke()
             <<"\n";
    }
#endif

    for(int i=s[ns].is(); i<=s[ns].ie(); i++) {
      for(int j=s[ns].js(); j<=s[ns].je(); j++) {
        for(int k=s[ns].ks(); k<=s[ns].ke(); k++) {
          if(clr->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
            real xcent=s[ns].x();
            real ycent=s[ns].y();
            real zcent=s[ns].z();
    
            real cseed = stratified_sphere(i,j,k, xcent,ycent,zcent);
    
            vol += clr->dV(i,j,k) * (1.0-cseed);
            if(clr->domain()->ibody().off(i,j,k-1)) {
              area_base += clr->dSz(i,j,k) * (1.0-cseed);
            }
          }
        }
      }
    }
  }

  boil::cart.sum_real(&vol);
  boil::cart.sum_real(&area_base);

  s[ns].set_vol_bubble(vol);
  s[ns].set_area_base(area_base);

#if 0
  if (ns==349) {
    boil::oout<<"vol_area:ns="<<ns<<" vol_bubble= "<<s[ns].vol_bubble()<<"\n";
    //exit(0);
  }
#endif

#ifdef DEBUG
  if (ns==0) {
    boil::aout<<"nucleation_vol_area: vol= "<<vol<<" area= "<<area_base<<"\n";
    boil::aout<<"nucleation_vol_area: rseed= "<<rseed<<" vol= "
              << 0.5*4.0/3.0*boil::pi*pow(rseed,3.0)<<" area= "
              << boil::pi*rseed*rseed<<"\n";
    boil::aout<<"nucleation_vol_area: sum_sink_energy= "<<s[ns].sum_sink_energy()<<"\n";
    //exit(0);
  }
#endif

  return;
}

/***************************************************************************//**
* add volume and area of dummy nucleation sites to genuine site
*******************************************************************************/
void Nucleation::vol_area_dummy(const int ns, const int nsd) {
  /* calculate volume and area of dummy site */
  vol_area(dsites);

  real vol = sites[ns].vol_bubble();  // volume of father
  real area = sites[ns].area_base();  // area of father
#if 0
  std::cout<<"vol_area_dummy: "<<ns<<" "<<nsd<<" vol= "<<vol<<" vol_dummy "<<dsites[nsd].vol_bubble()<<"\n";
#endif
  vol += dsites[nsd].vol_bubble();    // add dummy to father
  area += dsites[nsd].area_base();     // add dummy to father
  sites[ns].set_vol_bubble(vol);      // set to father
  sites[ns].set_area_base(area);      // set to father
  
  //dsites[nsd].set_vol_bubble(sites[ns].vol_bubble());
  //dsites[nsd].set_area_base (sites[ns].area_base() );
}
