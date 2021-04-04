#include "nucleation.h"

/***************************************************************************//**
*  plant single nucleation site
*******************************************************************************/
void Nucleation::plant_site(const int ns, const bool seed_source) {
  real vol_seed=0.0;   // volume of seed
  real area_base=0.0;  // area of bubble-base
  int kadj=0;          // k adjacent to wall
  
  for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
    for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
      for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
        //if (i<vf->si() || vf->ei()<i ||
        //    j<vf->sj() || vf->ej()<j) continue;
        if(vf->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
          real xcent=sites[ns].x();
          real ycent=sites[ns].y();
          real zcent=sites[ns].z();

          real cseed = stratified_sphere(i,j,k,
                                         xcent,ycent,zcent);

          if(matter_sig==Sign::pos()) {
            (*vf)[i][j][k]=std::min((*vf)[i][j][k],cseed);
          } else {
            (*vf)[i][j][k]=1.0-std::min(1.0-(*vf)[i][j][k],cseed);
          }

          vol_seed += vf->dV(i,j,k) * (1.0-cseed);
          if(vf->domain()->ibody().off(i,j,k-1)) {
            area_base += vf->dSz(Sign::neg(),i,j,k) * (1.0-cseed);
            kadj=k;
          }
        }
      }
    }
  }

  //std::cout<<"area_base= "<<area_base<<" "<<vol_seed<<" "<<kadj<<"\n";
  //std::cout<<"is,ie= "<<sites[ns].is()<<" "<<sites[ns].ie()<<"\n";

  if(qsrc&&seed_source) {

    /* set qsrc */
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        int k=kadj;
        real xcent=sites[ns].x();
        real ycent=sites[ns].y();
        real zcent=sites[ns].z();

        real cseed = stratified_sphere(i,j,k,
                                       xcent,ycent,zcent);
        if(matter_sig==Sign::neg()) {
          cseed = 1.-cseed;
        }
        if (vf->domain()->ibody().off(i,j,k-1)) {
          if (area_base>0.0) {
            (*qsrc)[i][j][k-1] -= rhov * vol_seed * latent
                               / (area_base * time->dt() * vf->dzc(k-1)) 
                               * (1.0-cseed) * vf->dV(i,j,k-1);
          }
        }
      } /* loop per site cells */
    } /* loop per site cells */
  } /* is there qsrc */

  return;
}

/***************************************************************************//**
*  plant single dummy site
*******************************************************************************/
void Nucleation::plant_dummy_site(const int nsd) {

  for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
    for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
      for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {

        real xcent=dsites[nsd].x();
        real ycent=dsites[nsd].y();
        real zcent=dsites[nsd].z();

        real cseed = stratified_sphere(i,j,k,
                                       xcent,ycent,zcent);

        if(matter_sig==Sign::pos()) {
          (*vf)[i][j][k]=std::min((*vf)[i][j][k],cseed);
        } else {
          (*vf)[i][j][k]=1.0-std::min(1.0-(*vf)[i][j][k],cseed);
        }

      } /* loop per dummy site cells */
    }
  }

  return;
}

