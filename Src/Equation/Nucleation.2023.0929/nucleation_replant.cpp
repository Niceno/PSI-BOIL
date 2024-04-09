#include "nucleation.h"
using namespace std;

/***************************************************************************//**
*  plant nucleation site
*******************************************************************************/
void Nucleation::replant () {

  if (bzoning!=true) {
    zoning();
  }

  /*-----------------+
  |  check criteria  |
  +-----------------*/

  real t_current = time->current_time();

  for(int ns=0; ns<size(); ns++){

    bool bseed=false;              // initialize bseed: bseed=true, if replant.
    real tpr_seed = tpr_site(ns);  // set seed temperature
    real clr_seed = clr_site(ns);  // set color function at seed point

    bool bheight = false;          // initialize bheight: bheight=true, if zplant is satisfied.
    if(sites[ns].zplant()<0.0) {   // crued code
      bheight = true;
    } else {
      real zft=zftmin( Range<real>(sites[ns].x()-dxmin,sites[ns].x()+dxmin)
                     , Range<real>(sites[ns].y()-dxmin,sites[ns].y()+dxmin)
                     , Range<real>(-boil::exa,boil::exa));
      if(zft>sites[ns].zplant()){
        bheight=true;
      }
    }

    /*---------------------------------------------------------------+
    |  replant if (1) seed temp. is higher than activation temp.     |
    |  and (2) did not replant in the last time step                 |
    |  and (3) seed point is liquid  (time_seed affect cutneck)      |
    |  and (4) bottom of previous bubble is higher than zplant()     |
    |  and (5) avoid replant immediately after cutneck.              |
    +---------------------------------------------------------------*/
    if( tpr_seed > sites[ns].active_tpr() && 
        sites[ns].seed_prev()==false &&
        clr_seed > 0.5  &&
        bheight &&
        t_current > (sites[ns].time_cutneck() + period_cut_replant) ) {

      /* check answers from outside of class to allow replant */
      if ( sites[ns].allow_replant() ){
        sites[ns].set_time_seed( t_current );
        bseed=true;
        boil::oout<<"replant:ns,x,y= "<<t_current<<" "<<ns<<" "
                  <<sites[ns].x()<<" "<<sites[ns].y()<<"\n";
      }

      /* request outside to allow replant */
      sites[ns].set_req_replant( true );
      //std::cout<<"replant: Pattern1 "<<boil::cart.iam()<<"\n";

    } else {
      sites[ns].set_req_replant( false );
      //std::cout<<"replant: Pattern2 "<<boil::cart.iam()<<"\n";
    }

    /*---------------------------------------+ 
    |  continue replant during seed_period   |
    +---------------------------------------*/
    if( t_current < (sites[ns].time_seed() + seed_period)) {
      bseed=true;
      sites[ns].set_req_replant( true );
      //std::cout<<"replant: Pattern2 "<<boil::cart.iam()<<"\n";
    }

    sites[ns].set_seed(bseed);
    //sites[ns].set_seed_prev(bseed);

    if (size()==1) {
      boil::oout<<"replant:printAll "<<ns<<" "<<t_current<<" "<<bseed<<" "
        <<tpr_seed<<" "<<clr_seed<<" "<<sites[ns].time_seed()<<" "
        <<sites[ns].time_cutneck()+period_cut_replant<<" "
        <<sites[ns].allow_replant()<<" "<<sites[ns].req_replant()<<"\n"; 
    }

  }

  /*--------------------------+
  |  copy bseed to dummy      |
  |  copy seed_prev to dummy  |
  +--------------------------*/
  for(int nsd=0; nsd<dsize(); nsd++){
    int fID=dsites[nsd].father();
    bool bseed = sites[fID].seed();
    dsites[nsd].set_seed(bseed);
    bool bseedp = sites[fID].seed_prev();
    dsites[nsd].set_seed_prev(bseedp);
  }

  /*----------+
  |  replant  |
  +----------*/
  real eps = dxmin * 1.5;
  real rhov = fluid()->rho(0);

  /* genuine sites */
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if (sites[ns].seed() ) {
      sites[ns].set_active(true);

      /* set color function */
      real vol_seed=0.0;   // volome of seed
      real area_base=0.0;  // area of bubble-base
      int kadj=0;          // k at adjasent to wall
      for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
        for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
          for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
              real xcent=sites[ns].x();
              real ycent=sites[ns].y();
              real zcent=sites[ns].z();
              real adist = sqrt( pow(clr->xc(i)-xcent,2.0)
                               + pow(clr->yc(j)-ycent,2.0)
                               + pow(clr->zc(k)-zcent,2.0))
                         - rseed;
              real cseed = 1.0;
              if(adist<-eps){
                cseed=0.0;
              }else if(adist<eps){
#if 0
                cseed= 0.5 + adist/(2.0*eps) 
                     + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
#else
                real mm=8;
                real x0=clr->xn(i);
                real y0=clr->yn(j);
                real z0=clr->zn(k);
                real ddx=clr->dxc(i)/real(mm);
                real ddy=clr->dyc(j)/real(mm);
                real ddz=clr->dzc(k)/real(mm);
                int itmp=0;
                for (int ii=0; ii<mm; ii++){
                  for (int jj=0; jj<mm; jj++){
                    for (int kk=0; kk<mm; kk++){
                      real xxc=x0+0.5*ddx+real(ii)*ddx;
                      real yyc=y0+0.5*ddy+real(jj)*ddy;
                      real zzc=z0+0.5*ddz+real(kk)*ddz;
                      real dist=sqrt(pow(xxc-xcent,2.0)
                                    +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
                      if (dist>rseed){
                        itmp=itmp+1;
                      }
                    }
                  }
                }
                cseed=real(itmp)/real(mm*mm*mm);
#endif
              }
              //(*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
              real clr_old = min(1.0,max(0.0,(*clr)[i][j][k]));
              real clr_new = std::min(clr_old,cseed);
              if(fabs(clr_new-0.5)<eps_clr) {
                clr_new = clrsurf+copysign(1.0,clr_new-clrsurf)*eps_clr;
              }
              (*clr)[i][j][k] = clr_new;
              vol_seed += ((1.0-clr_new)-(1.0-clr_old))*clr->dV(i,j,k);
              //vol_seed += clr->dV(i,j,k) * (1.0-cseed);

              if (clr->domain()->ibody().off(i,j,k-1)) {
                area_base += clr->dSz(i,j,k) * (1.0-cseed);
                kadj=k;
		// modify dmicro
                if( approx(sites[ns].time_seed(),t_current,boil::pico) ){
                  if(clr_new-clr_old<0.0){  // decrease-of-liquid (=increase-of-vapor)
                    dmicro[i][j][k]=min(dmicro0(i,j,k),dmicro[i][j][k]);  //update on 2021.12.18
                  }
                }
              }
            }
          }
        }
      }
#if 1
      /* set heat sink */
      if(heat_sink){
      if(sites[ns].seed_prev()==false) {
        boil::oout<<"##### nucl:replant USE HEAT-SINK "<<ns<<" #####"<<"\n";
        for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
          for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
            int k=kadj;
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
              real xcent=sites[ns].x();
              real ycent=sites[ns].y();
              real zcent=sites[ns].z();
              real adist = sqrt( pow(clr->xc(i)-xcent,2.0)
                               + pow(clr->yc(j)-ycent,2.0)
                               + pow(clr->zc(k)-zcent,2.0))
                         - rseed;
              real cseed = 1.0;
              if(adist<-eps){
                cseed=0.0;
              }else if(adist<eps){
#if 0
                cseed= 0.5 + adist/(2.0*eps) 
                     + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
#else
                real mm=8;
                real x0=clr->xn(i);
                real y0=clr->yn(j);
                real z0=clr->zn(k);
                real ddx=clr->dxc(i)/real(mm);
                real ddy=clr->dyc(j)/real(mm);
                real ddz=clr->dzc(k)/real(mm);
                int itmp=0;
                for (int ii=0; ii<mm; ii++){
                  for (int jj=0; jj<mm; jj++){
                    for (int kk=0; kk<mm; kk++){
                      real xxc=x0+0.5*ddx+real(ii)*ddx;
                      real yyc=y0+0.5*ddy+real(jj)*ddy;
                      real zzc=z0+0.5*ddz+real(kk)*ddz;
                      real dist=sqrt(pow(xxc-xcent,2.0)
                                    +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
                      if (dist>rseed){
                        itmp=itmp+1;
                      }
                    }
                  }
                }
                cseed=real(itmp)/real(mm*mm*mm);
#endif
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
      }
      }
#endif
    }
    sites[ns].set_seed_prev(sites[ns].seed());
  }

  /* dummy sites */
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    if (dsites[nsd].seed() ) {
      dsites[nsd].set_active(true);

      /* set color function */
      real vol_seed=0.0;   // volome of seed
      real area_base=0.0;  // area of bubble-base
      int kadj=0;          // k at adjasent to wall
      for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
        for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
          for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
              real xcent=dsites[nsd].x();
              real ycent=dsites[nsd].y();
              real zcent=dsites[nsd].z();
              real adist = sqrt( pow(clr->xc(i)-xcent,2.0)
                               + pow(clr->yc(j)-ycent,2.0)
                               + pow(clr->zc(k)-zcent,2.0))
                         - rseed;
              real cseed = 1.0;
              if(adist<-eps){
                cseed=0.0;
              }else if(adist<eps){
#if 0
      		      cseed= 0.5 + adist/(2.0*eps) 
                     + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
#else
                real mm=8;
                real x0=clr->xn(i);
                real y0=clr->yn(j);
                real z0=clr->zn(k);
                real ddx=clr->dxc(i)/real(mm);
                real ddy=clr->dyc(j)/real(mm);
                real ddz=clr->dzc(k)/real(mm);
                int itmp=0;
                for (int ii=0; ii<mm; ii++){
                  for (int jj=0; jj<mm; jj++){
                    for (int kk=0; kk<mm; kk++){
                      real xxc=x0+0.5*ddx+real(ii)*ddx;
                      real yyc=y0+0.5*ddy+real(jj)*ddy;
                      real zzc=z0+0.5*ddz+real(kk)*ddz;
                      real dist=sqrt(pow(xxc-xcent,2.0)
                                    +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
                      if (dist>rseed){
                        itmp=itmp+1;
                      }
                    }
                  }
                }
                cseed=real(itmp)/real(mm*mm*mm);
#endif
              }
              //(*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
              real clr_old = min(1.0,max(0.0,(*clr)[i][j][k]));
              real clr_new = std::min(clr_old,cseed);
              if(fabs(clr_new-0.5)<eps_clr) {
                clr_new = clrsurf+copysign(1.0,clr_new-clrsurf)*eps_clr;
              }
              (*clr)[i][j][k] = clr_new;
              vol_seed += ((1.0-clr_new)-(1.0-clr_old))*clr->dV(i,j,k);

              if (clr->domain()->ibody().off(i,j,k-1)) {
                area_base += clr->dSz(i,j,k) * (1.0-cseed);
                kadj=k;
                // modify dmicro
                int fID = dsites[nsd].father();
                if( approx(sites[fID].time_seed(),t_current,boil::pico) ){
                  if(clr_new-clr_old<0.0){  // decrease-of-liquid (=increase-of-vapor)
                    dmicro[i][j][k]=min(dmicro0(i,j,k),dmicro[i][j][k]);  //update on 2021.12.18
                  }
                }
              }
            }
          }
        }
      }
#if 1
      /* set heat sink */
      if(heat_sink){
      if(dsites[nsd].seed_prev()==false) {
        boil::oout<<"##### nucl:replant USE HEAT-SINK (dummy) #####"<<"\n";
        for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
          for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
            int k=kadj;
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
              real xcent=dsites[nsd].x();
              real ycent=dsites[nsd].y();
              real zcent=dsites[nsd].z();
              real adist = sqrt( pow(clr->xc(i)-xcent,2.0)
                               + pow(clr->yc(j)-ycent,2.0)
                               + pow(clr->zc(k)-zcent,2.0))
                         - rseed;
              real cseed = 1.0;
              if(adist<-eps){
                cseed=0.0;
              }else if(adist<eps){
#if 0
                cseed= 0.5 + adist/(2.0*eps)
                     + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
#else
                real mm=8;
                real x0=clr->xn(i);
                real y0=clr->yn(j);
                real z0=clr->zn(k);
                real ddx=clr->dxc(i)/real(mm);
                real ddy=clr->dyc(j)/real(mm);
                real ddz=clr->dzc(k)/real(mm);
                int itmp=0;
                for (int ii=0; ii<mm; ii++){
                  for (int jj=0; jj<mm; jj++){
                    for (int kk=0; kk<mm; kk++){
                      real xxc=x0+0.5*ddx+real(ii)*ddx;
                      real yyc=y0+0.5*ddy+real(jj)*ddy;
                      real zzc=z0+0.5*ddz+real(kk)*ddz;
                      real dist=sqrt(pow(xxc-xcent,2.0)
                                    +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
                      if (dist>rseed){
                        itmp=itmp+1;
                      }
                    }
                  }
                }
                cseed=real(itmp)/real(mm*mm*mm);
#endif
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
      }
      }
#endif
    }
    int fID = dsites[nsd].father();
    dsites[nsd].set_seed_prev(sites[fID].seed());
  }
  clr->exchange_all();

  st_active();

  /* set dmicro for new nucleation site */
#if 0
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if( approx(sites[ns].time_seed(),t_current,boil::pico) ){
      for(int k=sites[ns].ks(); k<=sites[ns].ke(); k++) {
        if (approx(clr->zn(k),zbtm,boil::pico)) {
          //cout<<boil::cart.iam()<<" "<<sites[ns].is()<<" "<<sites[ns].ie()<<" "<<sites[ns].js()<<" "<<sites[ns].je()<<" "<<clr->zn(k)<<"\n";
          for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
            for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
              if (i<clr->si() || clr->ei()<i ||
                  j<clr->sj() || clr->ej()<j) continue;
              if (area_vapor(i,j,k,Dir::kmin()) >0.0 ) {
                //dmicro[i][j][k]=dmicro0(i,j,k);
                dmicro[i][j][k]=min(dmicro0(i,j,k),dmicro[i][j][k]);  //update on 2021.12.17
                //cout<<i<<" "<<j<<" "<<k<<" "<<dmicro[i][j][k]<<"\n";
              }
            }
          }
        }
      }
    }
  }
#endif

  /* set dmicro for new dummy nucleation site */
#if 0
  //for(int nsd=0; nsd<dsize(); nsd++){
  for (int idd=0; idd<idd_nearRegion.size(); idd++){
    int nsd=idd_nearRegion[idd];
    int fID = dsites[nsd].father();
    if( approx(sites[fID].time_seed(),t_current,boil::pico) ){
      //cout<<"replant:k= "<<dsites[nsd].ks()<<" "<<dsites[nsd].ke()<<"\n";
      for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
        if (approx(clr->zn(k),zbtm,boil::pico)) {
          //cout<<boil::cart.iam()<<" "<<dsites[nsd].is()<<" "<<dsites[nsd].ie()<<" "<<dsites[nsd].js()<<" "<<dsites[nsd].je()<<" "<<clr->zn(k)<<"\n";
          for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
            for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
              if (i<clr->si() || clr->ei()<i ||
                  j<clr->sj() || clr->ej()<j) continue;
              if (area_vapor(i,j,k,Dir::kmin()) >0.0 ) {
                //dmicro[i][j][k]=dmicro0(i,j,k);
                dmicro[i][j][k]=min(dmicro0(i,j,k),dmicro[i][j][k]);  //update on 2021.12.17
                //cout<<i<<" "<<j<<" "<<k<<" "<<dmicro[i][j][k]<<"\n";
              }
            }
          }
        }
      }
    }
  }
#endif
}

