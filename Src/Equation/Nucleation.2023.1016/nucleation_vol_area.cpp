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
* copy volume and area of nucleation sites from genuine to dummy
*******************************************************************************/
void Nucleation::vol_area_dummy(const int ns, const int nsd) {
  dsites[nsd].set_vol_bubble(sites[ns].vol_bubble());
  dsites[nsd].set_area_base (sites[ns].area_base() );
}
