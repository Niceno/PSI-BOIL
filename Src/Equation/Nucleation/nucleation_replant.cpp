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

    bool bseed=false;              // bseed=true, if replant.
    real tpr_seed = tpr_site(ns);  // seed temperature
    real clr_seed = clr_site(ns);  // color function at seed point

    bool bheight = false;          // bheight=true, if zplant is satisfied.
    if(sites[ns].zplant()<0.0) {   // crude code
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
    bool clr_seed_cond;
    if(sig==Sign::pos()) {
      clr_seed_cond = clr_seed > 0.5;
    } else {
      clr_seed_cond = clr_seed < 0.5;
    }
    if( tpr_seed > sites[ns].active_tpr() && 
        sites[ns].seed_prev()==false &&
        clr_seed_cond &&
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

  /*----------------------+
  |  copy bseed to dummy  |
  +----------------------*/
  for(int nsd=0; nsd<dsize(); nsd++){
    int ns=dsites[nsd].father();
    bool bseed = sites[ns].seed();
    dsites[nsd].set_seed(bseed);
  }

  /*----------+
  |  replant  |
  +----------*/
  real eps = dxmin * 1.5;
  real rhov;
  if(sig==Sign::pos()) {
    rhov = fluid()->rho(0);
  } else {
    rhov = fluid()->rho(1);
  }


  /* genuine sites */
  //for(int ns=0; ns<size(); ns++){
  for (int id=0; id<id_nearRegion.size(); id++){
    int ns=id_nearRegion[id];
    if (sites[ns].seed() ) {
      sites[ns].set_active(true);
      /* set color function */
      real vol_seed=0.0;   // volome of seed
      real area_base=0.0;  // area of bubble-base
      int kadj=0;          // k adjacent to wall
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
              if(sig==Sign::pos()) {
                (*clr)[i][j][k]=std::min((*clr)[i][j][k],cseed);
              } else {
                (*clr)[i][j][k]=1.0-std::min(1.0-(*clr)[i][j][k],cseed);
              }
              vol_seed += clr->dV(i,j,k) * (1.0-cseed);
              if (clr->domain()->ibody().off(i,j,k-1)) {
                area_base += clr->dSz(Sign::neg(),i,j,k) * (1.0-cseed);
                kadj=k;
                //if(boil::cart.iam()==2)std::cout<<"area_base="<<area_base
                //   <<" "<<cseed<<"\n";
              }
            }
          }
        }
      }

#if 0
      /* set heat sink */
      if(sites[ns].seed_prev()==false) {
        for(int i=sites[ns].is(); i<=sites[ns].ie(); i++) {
          for(int j=sites[ns].js(); j<=sites[ns].je(); j++) {
            int k=kadj;
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
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
                //if(boil::cart.iam()==2)std::cout<<"area_base="<<area_base<<"\n";
                if (area_base>0.0) {
                  (*qsrc)[i][j][k-1] -= rhov * vol_seed * flu->latent(i,j,k)
                             / (area_base * time->dt() * clr->dzc(k-1))
                             * (1.0-cseed) * clr->dV(i,j,k-1);
                }
                //boil::aout<<"replant:heat sink "<<rhov<<" "<<vol_seed<<" "
                //  <<flu->latent(i,j,k)<<" "<<area_base<<" "<<time->dt()<<" "
                //  <<clr->dzc(k-1)<<" "<<(1.0-cseed)<<" "<<clr->dV(i,j,k-1)<<" "
                //  <<(*qsrc)[i][j][k-1]<<"\n";
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
      for(int i=dsites[nsd].is(); i<=dsites[nsd].ie(); i++) {
        for(int j=dsites[nsd].js(); j<=dsites[nsd].je(); j++) {
          for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
            if (i<clr->si() || clr->ei()<i ||
                j<clr->sj() || clr->ej()<j) continue;
            if( clr->domain()->ibody().on(i,j,k) ) {
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
    }
  }
  clr->exchange_all();

  st_active();

  /* set dmicro for new nucleation site */
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
      //cout<<"replant:k= "<<dsites[nsd].ks()<<" "<<dsites[nsd].ke()<<"\n";
      for(int k=dsites[nsd].ks(); k<=dsites[nsd].ke(); k++) {
        if (approx(clr->zn(k),zbtm,boil::pico)) {
          //cout<<boil::cart.iam()<<" "<<dsites[nsd].is()<<" "<<dsites[nsd].ie()<<" "<<dsites[nsd].js()<<" "<<dsites[nsd].je()<<" "<<clr->zn(k)<<"\n";
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

}
