#include "microlayer.h"
#include "../header.h"

/***************************************************************************//**
*  set dmicro for nucleation site
*******************************************************************************/
void Microlayer::upkeep_after_seeding() {

  std::cout<<"upkeep_after_seeding:begin\n";
  real t_current = time->current_time();

  /* set microlayer for genuine nucleation site */
  //for(int ns=0; ns<size(); ns++){
  for(int id=0; id<id_nearRegion.size(); id++) {
    int ns=id_nearRegion[id];
    //boil::oout<<"ML:upkeep= "<<sites[ns].time_seed()<<" "<<t_current<<"\n";
    if( approx(sites[ns].time_seed(),t_current,boil::pico) ){
      for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
        if (approx(vf->zn(k),zbtm,boil::pico)) {
          for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
            for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
              if (i<vf->si() || vf->ei()<i ||
                  j<vf->sj() || vf->ej()<j) continue;
              boil::oout<<"HERE "<<i<<" "<<k<<" "<<dmicro.dSz(Sign::neg(),i,j,k)<<" "<<area_vapor(Sign::neg(),Comp::k(),i,j,k)<<boil::endl;
#ifdef USE_VOF_NUCL
              if(in_vapor(i,j,k)) {
#else
              if(area_vapor(Sign::neg(),Comp::k(),i,j,k)>0.0) {
#endif
                //dmicro[i][j][k]=d0(i,j,k);
		dmicro[i][j][k]=min(d0(i,j,k),dmicro[i][j][k]); //Yohei 2022.07.24
              } else {
                dmicro[i][j][k]=boil::unreal;
              }
            }
          }
        }
      }
    }
  }

  /* set microlayer for dummy nucleation site */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    int fID = dsites[nsd].father();
    if( approx(sites[fID].time_seed(),t_current,boil::pico) ){
      for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
        if (approx(vf->zn(k),zbtm,boil::pico)) {
          for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
            for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
              if (i<vf->si() || vf->ei()<i ||
                  j<vf->sj() || vf->ej()<j) continue;
#ifdef USE_VOF_NUCL
              if(in_vapor(i,j,k)) {
#else
              if(area_vapor(Sign::neg(),Comp::k(),i,j,k)>0.0) {
#endif
                //dmicro[i][j][k]=d0(i,j,k);
		dmicro[i][j][k]=min(d0(i,j,k),dmicro[i][j][k]); //Yohei 2022.07.24
              } else {
                dmicro[i][j][k]=boil::unreal;
              }
            }
          }
        }
      }
    }
  }

  return;
}
