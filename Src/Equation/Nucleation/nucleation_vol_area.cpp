#include "nucleation.h"
//#define DEBUG

/***************************************************************************//**
* calculate volume and bottom_area of seed bubble
*******************************************************************************/
void Nucleation::vol_area(std::vector<Site> & s) {

  int ns=s.size()-1;

  real vol=0.0;   // volume of seed
  real area_base=0.0;  // area of bubble-base
  int kadj = 0;   // k in adjacent cell of wall (fluid side)
  
  if (sites[ns].contain_range() ) {
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
          //if (i<vf->si() || vf->ei()<i ||
          //    j<vf->sj() || vf->ej()<j) continue;
          if(vf->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
            real xcent=sites[ns].x();
            real ycent=sites[ns].y();
            real zcent=sites[ns].z();
    
            real cseed = stratified_sphere(i,j,k, xcent,ycent,zcent);
    
            vol += vf->dV(i,j,k) * (1.0-cseed);
            if(vf->domain()->ibody().off(i,j,k-1)) {
              area_base += vf->dSz(Sign::neg(),i,j,k) * (1.0-cseed);
	      kadj=k;
            }
          }
        }
      }
    }
  }

  boil::cart.sum_real(&vol);
  boil::cart.sum_real(&area_base);
  // kadj=0 if decomposed domain does not include near region around nucleation site

  s[ns].set_vol_bubble(vol);
  s[ns].set_area_base(area_base);
  s[ns].set_kadj(kadj);

#ifdef DEBUG
  if (ns==0) {
    boil::aout<<"nucleation_vol_area: vol= "<<vol<<" area= "<<area_base
              <<" kadj= "<<kadj<<"\n";
    boil::aout<<"nucleation_vol_area: rseed= "<<rseed<<" vol= "
              << 0.5*4.0/3.0*boil::pi*pow(rseed,3.0)<<" area= "
              << boil::pi*rseed*rseed<<"\n";
    boil::aout<<"nucleation_vol_area: sum_sink_energy= "<<s[ns].sum_sink_energy()<<"\n";
    //exit(0);
  }
#endif

  return;
}
