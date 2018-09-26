#include "nucleation.h"

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::plant () {

  real rhov = fluid()->rho(0);

  if (bzoning!=true) {
    zoning();
  }

  real eps = dxmin * 1.5;
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];

    //boil::oout<<"nucl.plant: ID"<<ns<<" at: "<<sites[ns].x()<<", "
    //          <<sites[ns].y()<<", "<<sites[ns].z()<<"\n";

    sites[ns].set_time_seed( time->current_time() );
    sites[ns].set_active(true);
    /* set color function */
    real vol_seed=0.0;   // volome of seed
    real area_base=0.0;  // area of bubble-base
    int kadj=0;          // k at adjasent to wall
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
          if(clr->domain()->ibody().fV(i,j,k)!=0.0) {  // R4
            real adist = sqrt( pow(clr->xc(i)-sites[ns].x(),2.0)
                             + pow(clr->yc(j)-sites[ns].y(),2.0)
                             + pow(clr->zc(k)-sites[ns].z(),2.0))
                       - rseed;
            real cseed = 1.0;
            if(adist<-eps){
              cseed=0.0;
            }else if(adist<eps){
              cseed= 0.5 + adist/(2.0*eps) 
                   + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
            }
            (*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
            vol_seed += clr->dV(i,j,k) * (1.0-cseed);
            if (clr->domain()->ibody().off(i,j,k-1)) {
              area_base += clr->dSz(i,j,k) * (1.0-cseed);
              kadj=k;
            }
          }
        }
      }
    }
    //std::cout<<"area_base= "<<area_base<<" "<<vol_seed<<" "<<kadj<<"\n";
    //std::cout<<"is,ie= "<<sites[ns].is()<<" "<<sites[ns].ie()<<"\n";

    /* set qsrc */
    for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
      for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
        int k=kadj;
        real adist = sqrt( pow(clr->xc(i)-sites[ns].x(),2.0)
                         + pow(clr->yc(j)-sites[ns].y(),2.0)
                         + pow(clr->zc(k)-sites[ns].z(),2.0))
                   - rseed;
        real cseed = 1.0;
        if(adist<-eps){
          cseed=0.0;
        }else if(adist<eps){
          cseed= 0.5 + adist/(2.0*eps) 
               + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
        }
        if (clr->domain()->ibody().off(i,j,k-1)) {
          if (area_base>0.0) {
            (*qsrc)[i][j][k-1] -= rhov * vol_seed * latent
                               / (area_base * time->dt() * clr->dzc(k-1)) 
                               * (1.0-cseed) * clr->dV(i,j,k-1);
          }
        }
      }
    }
  }
  qsrc->exchange();

  /* dummy sites */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    //boil::oout<<"nucl.plantdummy: ID"<<nsd<<" at: "<<dsites[nsd].x()<<", "
    //          <<dsites[nsd].y()<<", "<<dsites[nsd].z()<<"\n";
    dsites[nsd].set_active(true);
    for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
      for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
        for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
          real adist = sqrt( pow(clr->xc(i)-dsites[nsd].x(),2.0)
                           + pow(clr->yc(j)-dsites[nsd].y(),2.0)
                           + pow(clr->zc(k)-dsites[nsd].z(),2.0))
                     - rseed;
          real cseed = 1.0;
          if(adist<-eps){
            cseed=0.0;
          }else if(adist<eps){
            cseed= 0.5 + adist/(2.0*eps) 
                 + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
          }
          (*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
        }
      }
    }
  }

  clr->exchange_all();

  st_active();

  real t_current = time->current_time();

  /* set dmicro for new nucleation site */
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    //std::cout<<sites[ns].time_seed()<<" "<<t_current<<"\n";
    if( approx(sites[ns].time_seed(),t_current,boil::pico) ){
      for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
        if (approx(clr->zn(k),zbtm,boil::pico)) {
          for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
            for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
              if (i<clr->si() || clr->ei()<i ||
                  j<clr->sj() || clr->ej()<j) continue;
              if (area_vapor(i,j,k,Dir::kmin()) >0.0 ) {
                dmicro[i][j][k]=dmicro0(i,j,k);
                //cout<<i<<" "<<j<<" "<<k<<" "<<dmicro[i][j][k]<<"\n";
              }
            }
          }
        }
      }
    }
  }

  /* set dmicro for new dummy nucleation site */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    int fID = dsites[nsd].father();
    if( approx(sites[fID].time_seed(),t_current,boil::pico) ){
      for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
        if (approx(clr->zn(k),zbtm,boil::pico)) {
          for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
            for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
              if (i<clr->si() || clr->ei()<i ||
                  j<clr->sj() || clr->ej()<j) continue;
              if (area_vapor(i,j,k,Dir::kmin()) >0.0 ) {
                dmicro[i][j][k]=dmicro0(i,j,k);
                //cout<<i<<" "<<j<<" "<<k<<" "<<dmicro[i][j][k]<<"\n";
              }
            }
          }
        }
      }
    }
  }

#if 0
  conc.insert_bc(c);
  c.exchange_all();
#endif
}

/*-----------------------------------------------------------------------------+
 '$Id: nucleation_plant.cpp,v 1.7 2018/09/17 12:47:36 sato Exp $'/
+-----------------------------------------------------------------------------*/
