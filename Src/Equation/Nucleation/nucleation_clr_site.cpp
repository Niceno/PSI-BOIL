#include "nucleation.h"

/******************************************************************************/
real Nucleation::clr_site(const int ns ) {
/***************************************************************************//**
*  \brief return color function just above the site
*******************************************************************************/

  real clr_seed = -boil::unreal;
  if (sites[ns].contain_site()) {
    clr_seed = (*clr)[sites[ns].ic_site()]
                     [sites[ns].jc_site()]
                     [sites[ns].kc_site()+1];  // fluid phase, crude code
  }
  boil::cart.max_real(&clr_seed);

  if(!boil::realistic(clr_seed)) {
    boil::oout<<"nucleation:clr_site: Error!!! Cannot find seed point!\n";
    exit(0);
  }

  // if the clr just above nucleation site is vapor, then stop investigation
  //if (in_vapor(clr_seed)) {
  if (clr_seed < threshold_clr_site) {
    return clr_seed;
  }

  // clr just above nucleation site is in liquid
  // check the minimu and max of the region where bubble will be plant
  real clr_min =  boil::unreal;
  real clr_max = -boil::unreal;
  real dxy = 0.0;  // grid spacing at nucleation site
  if(sites[ns].contain_site()){
    real dxx = (*clr).dxc(sites[ns].ic_site());
    real dyy = (*clr).dyc(sites[ns].jc_site());
    dxy = boil::maxr(dxx,dyy);
  }
  boil::cart.max_real(&dxy);
  if (dxy==0.0) {
    boil::oout<<"nucleation_clr_site:Error dxy=0.0\n";
    exit(0);
  }

  if (sites[ns].contain_range()) {
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
          if(clr->domain()->ibody().fV(i,j,k)!=0.0) {
            real xcent = sites[ns].x();
            real ycent = sites[ns].y();
            real zcent = sites[ns].z();
	    real xx = (*clr).xc(i);
	    real yy = (*clr).yc(j);
	    real zz = (*clr).zc(k);
            real dist = sqrt(pow(xx-xcent,2.0)
			    +pow(yy-ycent,2.0)
			    +pow(zz-zcent,2.0));
            //if (dist < rseed) {
            if (dist < rseed + rseed_plus*dxy) {  // rseed_plus cell outside
              if (clr_min>(*clr)[i][j][k]) clr_min=(*clr)[i][j][k];
              if (clr_max<(*clr)[i][j][k]) clr_max=(*clr)[i][j][k];
	    }
          }
        }
      }
    }
  }
  boil::cart.min_real(&clr_min);
  boil::cart.max_real(&clr_max);

  return clr_min;
}
