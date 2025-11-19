#include "body.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"
#include "../Domain/domain.h"
#include <iomanip>
//#define DEBUG

/******************************************************************************/
void Body::distfunc(const Domain & dom) {
/***************************************************************************//**
*  /brief calculate signed distance function from sca.
*         Reference:G.,Russo and P.,Smereka,"A remark on computing distance
*         functions",J.Comp.Phys.,Vol.163,2000,pp.51-67
*           input    :sca
*           output   :bdist
*           temporary:stmp,dflag
*******************************************************************************/
#ifdef DEBUG
  std::cout<<"body_distfunc::start "<<boil::cart.iam()<<"\n";
#endif

  /*------------------------------+
  |  calculate distance function  |
  +------------------------------*/
  int irank=boil::cart.iam();
  int nproc=boil::cart.nproc();
  int ipcal[nproc];
  for(int i=0; i<nproc; i++){
    ipcal[i]=0;
  }
  /* define scalar */
  stmp   = new Scalar(dom); (*stmp ) = 0.0; (*stmp )=sca->shape();
  dflag  = new Scalar(dom); (*dflag) = 0.0; (*dflag)=sca->shape();
  for_avijk((*bdist),i,j,k){
    (*bdist)[i][j][k]=1.0e+60;
  }

#ifdef DEBUG
  std::cout<<"distfunc::nccells= "<<nccells()<<" "<<irank<<"\n";
  std::cout<<"distfunc::ccells.size = "<<cells[3].size()<<" "<<irank<<"\n";
#endif
  for(int cc=0; cc<cells[3].size(); cc++) {
    int ic, jc, kc;
    dom.ibody().ijk(cc,&ic,&jc,&kc);
    for(int i=-1; i<=1; i++){
      for(int j=-1; j<=1; j++){
        for(int k=-1; k<=1; k++){
          if( ic+i<0 || ic+i>=bdist->ni()) break;
          if( jc+j<0 || jc+j>=bdist->nj()) break;
          if( kc+k<0 || kc+k>=bdist->nk()) break;
          real dd;
          for(int iv=0; iv<polys[cc].nnodes(); iv++){
            dd = sqrt( pow( polys[cc].xn(iv) - bdist->xc(ic+i) , 2.0)
                     + pow( polys[cc].yn(iv) - bdist->yc(jc+j) , 2.0)
                     + pow( polys[cc].zn(iv) - bdist->zc(kc+k) , 2.0));
            if(dd<fabs((*bdist)[ic+i][jc+j][kc+k])){
              if(fV(ic+i,jc+j,kc+k)>0.5){
                (*bdist)[ic+i][jc+j][kc+k] =  dd;
              } else {
                (*bdist)[ic+i][jc+j][kc+k] = -dd;
              }
            }
          }
          /* parpendicular line */
          if(polys[cc].foot( bdist->xc(ic+i)
                           , bdist->yc(jc+j), bdist->zc(kc+k),dd)){
            if(dd<fabs((*bdist)[ic+i][jc+j][kc+k])){
              if(fV(ic+i,jc+j,kc+k)>0.5){
                (*bdist)[ic+i][jc+j][kc+k] =  dd;
              } else {
                (*bdist)[ic+i][jc+j][kc+k] = -dd;
              }
            }
          }
        }
      }
    }
  }

  if(cells[3].size()>=1)ipcal[irank]=1;
  boil::cart.sum_int_n(ipcal,nproc);
  insert_bc_dist(bdist);
  bdist->exchange(ipcal);

  /* set constants */
  int nlayer = 12;
  const real distmin = -dxmin * nlayer;
  const real distmax =  dxmin * nlayer;
  const int nflagmax =  nlayer+4;
  const int nflagmin = -nlayer-4;

  /* set dflag */
  for_avijk((*dflag),i,j,k) {
    if((*sca)[i][j][k]<0.5){
      (*dflag)[i][j][k]=nflagmin;
    } else {
      (*dflag)[i][j][k]=nflagmax;
    }
  }

  /* including or next to free-surface (NFCell) */
  /* i-direction */
  for(int i=sifl[3]; i<eifl[3]; i++){
    for(int j=sjfl[3]; j<=ejfl[3]; j++){
      for(int k=skfl[3]; k<=ekfl[3]; k++){
        if(((*sca)[i][j][k]-0.5)*((*sca)[i+1][j][k]-0.5)<=0.0){
           (*dflag)[i  ][j][k]=0;
           (*dflag)[i+1][j][k]=0;
        }
      }
    }
  }
  /* j-direction */
  for(int i=sifl[3]; i<=eifl[3]; i++){
    for(int j=sjfl[3]; j<ejfl[3]; j++){
      for(int k=skfl[3]; k<=ekfl[3]; k++){
        if(((*sca)[i][j][k]-0.5)*((*sca)[i][j+1][k]-0.5)<=0.0){
          (*dflag)[i][j  ][k]=0;
          (*dflag)[i][j+1][k]=0;
        }
      }
    }
  }
  /* k-direction */
  for(int i=sifl[3]; i<=eifl[3]; i++){
    for(int j=sjfl[3]; j<=ejfl[3]; j++){
      for(int k=skfl[3]; k<ekfl[3]; k++){
        if(((*sca)[i][j][k]-0.5)*((*sca)[i][j][k+1]-0.5)<=0.0){
          (*dflag)[i][j][k  ]=0;
          (*dflag)[i][j][k+1]=0;
        }
      }
    }
  }
  dflag->exchange(ipcal);

  int idelay;
  if(ipcal[irank]==1)idelay=0;
  if(ipcal[irank]==0)idelay=0;

  for(int layer=1; layer<=nlayer; layer++){
    /* i-direction */
    for(int i=sifl[3]; i<eifl[3]; i++){
      for(int j=sjfl[3]; j<=ejfl[3]; j++){
        for(int k=skfl[3]; k<=ekfl[3]; k++){
          if(fabs((*dflag)[i  ][j][k])==nflagmax &&
             fabs((*dflag)[i+1][j][k])==(layer-1-idelay)) {
            (*dflag)[i  ][j][k] = layer*int(copysign(1.0,(*dflag)[i  ][j][k]));
            if((*bdist)[i][j][k]==0.0) (*bdist)[i][j][k]=(*bdist)[i+1][j][k];
          } else if(fabs((*dflag)[i+1][j][k])==nflagmax &&
                    fabs((*dflag)[i  ][j][k])==(layer-1-idelay)){
            (*dflag)[i+1][j][k] = layer*int(copysign(1.0,(*dflag)[i+1][j][k]));
            if((*bdist)[i+1][j][k]==0.0) (*bdist)[i+1][j][k]=(*bdist)[i][j][k];
          }
        }
      }
    }
    /* j-direction */
    for(int i=sifl[3]; i<=eifl[3]; i++){
      for(int j=sjfl[3]; j<ejfl[3]; j++){
        for(int k=skfl[3]; k<=ekfl[3]; k++){
          if(fabs((*dflag)[i][j  ][k])==nflagmax &&
             fabs((*dflag)[i][j+1][k])==(layer-1)) {
            (*dflag)[i  ][j][k] = layer*int(copysign(1.0,(*dflag)[i  ][j][k]));
          } else if(fabs((*dflag)[i][j+1][k])==nflagmax &&
                    fabs((*dflag)[i][j  ][k])==(layer-1)) {
            (*dflag)[i][j+1][k] = layer*int(copysign(1.0,(*dflag)[i][j+1][k]));
          }
        }
      }
    }
    /* k-direction */
    for(int i=sifl[3]; i<=eifl[3]; i++){
      for(int j=sjfl[3]; j<=ejfl[3]; j++){
        for(int k=skfl[3]; k<ekfl[3]; k++){
          if(fabs((*dflag)[i][j][k  ])==nflagmax &&
             fabs((*dflag)[i][j][k+1])==(layer-1)) {
            (*dflag)[i  ][j][k] = layer*int(copysign(1.0,(*dflag)[i  ][j][k]));
          } else if(fabs((*dflag)[i][j][k+1])==nflagmax &&
                    fabs((*dflag)[i][j][k  ])==(layer-1)) {
            (*dflag)[i][j][k+1] = layer*int(copysign(1.0,(*dflag)[i][j][k+1]));
          }
        }
      }
    }
    dflag->exchange(ipcal);
  }
  dflag->exchange();
  insert_bc_dist((dflag));

  /*-----------------------------------------+
  |  initial condition of distance function  |
  +-----------------------------------------*/
  for_avijk((*dflag),i,j,k){
    if((*dflag)[i][j][k]!=0){
      (*bdist)[i][j][k] = dxmin * (*dflag)[i][j][k];
    }
  }

  /*-----------------------+
  |  distance calculation  |
  +-----------------------*/
  for(int it=1; it<=24; it++){
    for_vijk((*bdist),i,j,k){
      if(fabs((*dflag)[i][j][k])<=nlayer){
        if((*dflag)[i][j][k]==0.0){
          (*stmp)[i][j][k]=0.0;
        } else {
          real dx = std::min(bdist->dzc(k),
                    std::min(bdist->dxc(i),bdist->dyc(j))); // crude
          real dtau = 0.5*dx;
          real asgn=copysign(1.0, (*bdist)[i][j][k]);
          /* cells except NFCell */
          real a=((*bdist)[i][j][k]-(*bdist)[i-1][j][k])/bdist->dxw(i);
          real b=((*bdist)[i+1][j][k]-(*bdist)[i][j][k])/bdist->dxe(i);
          real c=((*bdist)[i][j][k]-(*bdist)[i][j-1][k])/bdist->dys(j);
          real d=((*bdist)[i][j+1][k]-(*bdist)[i][j][k])/bdist->dyn(j);
          real e=((*bdist)[i][j][k]-(*bdist)[i][j][k-1])/bdist->dzb(k);
          real f=((*bdist)[i][j][k+1]-(*bdist)[i][j][k])/bdist->dzt(k);
          real g;
          if((*bdist)[i][j][k]>=0){
            real ap = std::max(a,0.0);
            real bm = std::min(0.0,b);
            real cp = std::max(c,0.0);
            real dm = std::min(0.0,d);
            real ep = std::max(e,0.0);
            real fm = std::min(0.0,f);
            g = sqrt(std::max(ap*ap,bm*bm)
                    +std::max(cp*cp,dm*dm)
                    +std::max(ep*ep,fm*fm))-1.0;
          } else {
            real am = std::min(0.0,a);
            real bp = std::max(b,0.0);
            real cm = std::min(0.0,c);
            real dp = std::max(d,0.0);
            real em = std::min(0.0,e);
            real fp = std::max(f,0.0);
            g = sqrt(std::max(am*am,bp*bp)
                    +std::max(cm*cm,dp*dp)
                    +std::max(em*em,fp*fp))-1.0;
          }
          (*stmp)[i][j][k] = dtau * asgn * g;
        }
      }
    }

      /* update Gauss-Saidel */
    for_vijk((*bdist),i,j,k){
      (*bdist)[i][j][k] -= (*stmp)[i][j][k];
    }

    insert_bc_dist(bdist);
    (*bdist).exchange_all();
  }

  /* release memory */
  delete stmp,dflag;

  //boil::plot->plot(*bdist,"bdist");

  return;
}
