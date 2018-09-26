#include "cipcsl2.h"
#include "../../../Parallel/Out/out.h"
#include <iomanip>
#define WALL
using namespace std;

/******************************************************************************/
void CIPCSL2::set_alp() {
/***************************************************************************//**
*  \brief set alp for redistance.
*         output : alp
*         tempolary: stmp
*******************************************************************************/

#ifdef DEBUG
  boil::oout<<"cipcsl2_set_alp:start\n";
#endif


if (!localSharpen) {
  /* global sharpening */
  if(ialpcal==0)
    for_aijk(i,j,k)
      alp[i][j][k]=0;
  ialpcal=1;
} else {
  /* local sharpening */
  /*-----------+
  |  set flag  |
  +-----------*/
  /* initialize */
  for_aijk(i,j,k)
    stmp[i][j][k]=0.0;
  for_aijk(i,j,k)
    alp[i][j][k]=0.0;


  /* i-direction */
  for(int i=si()-1; i<=ei(); i++)
    for_jk(j,k)
      if((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0){
        stmp[i  ][j][k]=1.0;
        stmp[i+1][j][k]=1.0;
      }

  /* j-direction */
  for(int j=sj()-1; j<=ej(); j++)
    for_ik(i,k)
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        stmp[i][j  ][k]=1.0;
        stmp[i][j+1][k]=1.0;
      }

  /* k-direction */
  for(int k=sk()-1; k<=ek(); k++)
    for_ij(i,j)
      if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        stmp[i][j][k  ]=1.0;
        stmp[i][j][k+1]=1.0;
      }

  /* immersed boundary */
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k){
      if(dom->ibody().off(i,j,k))
        stmp[i][j][k]=0.0;
    }
  }

  stmp.exchange();

  /*------------+
  |  cal index  |
  +------------*/
  real aindex_max=0.0;
  for_ijk(i,j,k){
    if(stmp[i][j][k]!=0.0){
      real dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz, c1,c2,c3;
      real aindex;
      dudx=((*u)[Comp::u()][i+1][j][k]-(*u)[Comp::u()][i][j][k])/clr.dxc(i);
      dudy=0.5*((*u)[Comp::u()][i  ][j+1][k]-(*u)[Comp::u()][i  ][j-1][k]
               +(*u)[Comp::u()][i+1][j+1][k]-(*u)[Comp::u()][i+1][j-1][k])
          / (clr.dys(j)+clr.dyn(j));
      dudz=0.5*((*u)[Comp::u()][i  ][j][k+1]-(*u)[Comp::u()][i  ][j][k-1]
               +(*u)[Comp::u()][i+1][j][k+1]-(*u)[Comp::u()][i+1][j][k-1])
          / (clr.dzb(k)+clr.dzt(k));
      dvdx=0.5*((*u)[Comp::v()][i+1][j  ][k]-(*u)[Comp::v()][i-1][j  ][k]
               +(*u)[Comp::v()][i+1][j+1][k]-(*u)[Comp::v()][i-1][j+1][k])
          / (clr.dxw(i)+clr.dxe(i));
      dvdy=((*u)[Comp::v()][i][j+1][k]-(*u)[Comp::v()][i][j][k])/clr.dyc(j);
      dvdz=0.5*((*u)[Comp::v()][i][j  ][k+1]-(*u)[Comp::v()][i][j  ][k-1]
               +(*u)[Comp::v()][i][j+1][k+1]-(*u)[Comp::v()][i][j+1][k-1])
          / (clr.dzb(k)+clr.dzt(k));
      dwdx=0.5*((*u)[Comp::w()][i+1][j][k  ]-(*u)[Comp::w()][i-1][j][k  ]
               +(*u)[Comp::w()][i+1][j][k+1]-(*u)[Comp::w()][i-1][j][k+1])
          / (clr.dxw(i)+clr.dxe(i));
      dwdy=0.5*((*u)[Comp::w()][i][j+1][k  ]-(*u)[Comp::w()][i][j-1][k  ]
               +(*u)[Comp::w()][i][j+1][k+1]-(*u)[Comp::w()][i][j-1][k+1])
          / (clr.dys(j)+clr.dyn(j));
      dwdz=((*u)[Comp::w()][i][j][k+1]-(*u)[Comp::w()][i][j][k])/clr.dzc(k);

      c1 = dudx*nx[i][j][k] + dvdx*ny[i][j][k] + dwdx*nz[i][j][k];
      c2 = dudy*nx[i][j][k] + dvdy*ny[i][j][k] + dwdy*nz[i][j][k];
      c3 = dudz*nx[i][j][k] + dvdz*ny[i][j][k] + dwdz*nz[i][j][k];
      aindex= (c1*nx[i][j][k]+c2*ny[i][j][k]+c3*nz[i][j][k]);
      alp[i][j][k]=aindex;
      aindex_max=std::max(aindex,aindex_max);
    }
  }

  /* cells next to immersed boundary */
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    if(stmp[i][j][k]==0.0) continue;

    real dudx,dudy,dudz, dvdx,dvdy,dvdz, dwdx,dwdy,dwdz, c1,c2,c3;
    real aindex;
    real dxw, dxe, dys, dyn, dzb, dzt;

    // u
    dxw = 0.5*clr.dxc(i);
    dxe = 0.5*clr.dxc(i);
    dys = clr.dys(j);
    dyn = clr.dyn(j);
    dzb = clr.dzb(k);
    dzt = clr.dzt(k);
    real uw = (*u)[Comp::u()][i  ][j][k];
    real ue = (*u)[Comp::u()][i+1][j][k];
    real us = 0.5*((*u)[Comp::u()][i  ][j-1][k]+(*u)[Comp::u()][i+1][j-1][k]);
    real un = 0.5*((*u)[Comp::u()][i  ][j+1][k]+(*u)[Comp::u()][i+1][j+1][k]);
    real ub = 0.5*((*u)[Comp::u()][i  ][j][k-1]+(*u)[Comp::u()][i+1][j][k-1]);
    real ut = 0.5*((*u)[Comp::u()][i  ][j][k+1]+(*u)[Comp::u()][i+1][j][k+1]);
    if(dom->ibody().off(Comp::u(),i,j,k)){
      dxw = clr.dxw(i)*dom->ibody().fdxw(i,j,k);
      uw = 0.0;
    }
    if(dom->ibody().off(Comp::u(),i+1,j,k)){
      dxe = clr.dxe(i)*dom->ibody().fdxe(i,j,k);
      ue = 0.0;
    }
    if(dom->ibody().fdys(i,j,k)<1.0){
      dys *= dom->ibody().fdys(i,j,k);
      us = 0.0;
    }
    if(dom->ibody().fdyn(i,j,k)<1.0){
      dyn *= dom->ibody().fdyn(i,j,k);
      un = 0.0;
    }
    if(dom->ibody().fdzb(i,j,k)<1.0){
      dzb *= dom->ibody().fdzb(i,j,k);
      ub = 0.0;
    }
    if(dom->ibody().fdzt(i,j,k)<1.0){
      dzt *= dom->ibody().fdzt(i,j,k);
      ut = 0.0;
    }
    dudx=(ue-uw)/(dxw+dxe);
    dudy=(un-us)/(dys+dyn);
    dudz=(ut-ub)/(dzb+dzt);

    // v
    dxw = clr.dxw(i);
    dxe = clr.dxe(i);
    dys = 0.5*clr.dyc(j);
    dyn = 0.5*clr.dyc(j);
    dzb = clr.dzb(k);
    dzt = clr.dzt(k);
    real vw = 0.5*((*u)[Comp::v()][i-1][j  ][k]+(*u)[Comp::v()][i-1][j+1][k]);
    real ve = 0.5*((*u)[Comp::v()][i+1][j  ][k]+(*u)[Comp::v()][i+1][j+1][k]);
    real vs = (*u)[Comp::v()][i][j  ][k];
    real vn = (*u)[Comp::v()][i][j+1][k];
    real vb = 0.5*((*u)[Comp::v()][i][j  ][k-1]+(*u)[Comp::v()][i][j+1][k-1]);
    real vt = 0.5*((*u)[Comp::v()][i][j  ][k+1]+(*u)[Comp::v()][i][j+1][k+1]);
    if(dom->ibody().fdxw(i,j,k)<1.0){
      dxw *= dom->ibody().fdxw(i,j,k);
      vw = 0.0;
    }
    if(dom->ibody().fdxe(i,j,k)<1.0){
      dxe *= dom->ibody().fdxe(i,j,k);
      ve = 0.0;
    }
    if(dom->ibody().off(Comp::v(),i,j,k)){
      dys = clr.dys(j)*dom->ibody().fdys(i,j,k);
      vs = 0.0;
    }
    if(dom->ibody().off(Comp::v(),i,j+1,k)){
      dyn = clr.dyn(j)*dom->ibody().fdyn(i,j,k);
      vn = 0.0;
    }
    if(dom->ibody().fdzb(i,j,k)<1.0){
      dzb *= dom->ibody().fdzb(i,j,k);
      vb = 0.0;
    }
    if(dom->ibody().fdzt(i,j,k)<1.0){
      dzt *= dom->ibody().fdzt(i,j,k);
      vt = 0.0;
    }
    dvdx=0.5*((*u)[Comp::v()][i+1][j  ][k]-(*u)[Comp::v()][i-1][j  ][k]
             +(*u)[Comp::v()][i+1][j+1][k]-(*u)[Comp::v()][i-1][j+1][k])
        / (clr.dxw(i)+clr.dxe(i));
    dvdy=((*u)[Comp::v()][i][j+1][k]-(*u)[Comp::v()][i][j][k])/clr.dyc(j);
    dvdz=0.5*((*u)[Comp::v()][i][j  ][k+1]-(*u)[Comp::v()][i][j  ][k-1]
             +(*u)[Comp::v()][i][j+1][k+1]-(*u)[Comp::v()][i][j+1][k-1])
        / (clr.dzb(k)+clr.dzt(k));
    dvdx=(ve-vw)/(dxw+dxe);
    dvdy=(vn-vs)/(dys+dyn);
    dvdz=(vt-vb)/(dzb+dzt);

    // w
    dxw = clr.dxw(i);
    dxe = clr.dxe(i);
    dys = clr.dys(j);
    dyn = clr.dyn(j);
    dzb = 0.5*clr.dzc(k);
    dzt = 0.5*clr.dzc(k);
    real ww = 0.5*((*u)[Comp::w()][i-1][j][k  ]+(*u)[Comp::w()][i-1][j][k+1]);
    real we = 0.5*((*u)[Comp::w()][i+1][j][k  ]+(*u)[Comp::w()][i+1][j][k+1]);
    real ws = 0.5*((*u)[Comp::w()][i][j-1][k  ]+(*u)[Comp::w()][i][j-1][k+1]);
    real wn = 0.5*((*u)[Comp::w()][i][j+1][k  ]+(*u)[Comp::w()][i][j+1][k+1]);
    real wb = (*u)[Comp::w()][i][j][k];
    real wt = (*u)[Comp::w()][i][j][k+1];
    if(dom->ibody().fdxw(i,j,k)<1.0){
      dxw *= dom->ibody().fdxw(i,j,k);
      ww = 0.0;
    }
    if(dom->ibody().fdxe(i,j,k)<1.0){
      dxe *= dom->ibody().fdxe(i,j,k);
      we = 0.0;
    }
    if(dom->ibody().fdys(i,j,k)<1.0){
      dys *= dom->ibody().fdys(i,j,k);
      ws = 0.0;
    }
    if(dom->ibody().fdyn(i,j,k)<1.0){
      dyn *= dom->ibody().fdyn(i,j,k);
      wn = 0.0;
    }
    if(dom->ibody().off(Comp::w(),i,j,k)){
      dzb = clr.dzb(k)*dom->ibody().fdzb(i,j,k);
      wb = 0.0;
    }
    if(dom->ibody().off(Comp::w(),i,j,k+1)){
      dzt = clr.dzt(k)*dom->ibody().fdzt(i,j,k);
      wt = 0.0;
    }
    dwdx=(we-ww)/(dxw+dxe);
    dwdy=(wn-ws)/(dys+dyn);
    dwdz=(wt-wb)/(dzb+dzt);

    c1 = dudx*nx[i][j][k] + dvdx*ny[i][j][k] + dwdx*nz[i][j][k];
    c2 = dudy*nx[i][j][k] + dvdy*ny[i][j][k] + dwdy*nz[i][j][k];
    c3 = dudz*nx[i][j][k] + dvdz*ny[i][j][k] + dwdz*nz[i][j][k];
    aindex= (c1*nx[i][j][k]+c2*ny[i][j][k]+c3*nz[i][j][k]);
    alp[i][j][k]=aindex;
    aindex_max=std::max(aindex,aindex_max);
  }

  insert_bc_alp(alp);
  ib_ext_scalar(alp);
  alp.exchange();

  /*--------------+
  |  average alp  |
  +--------------*/
  /* i-direction */
  fn=0.0;
  for(int i=si()-1; i<=ei(); i++)
    for_jk(j,k)
      if((clr[i][j][k]-phisurf)*(clr[i+1][j][k]-phisurf)<=0.0){
        atmp[i  ][j][k]=alp[i][j][k]+alp[i+1][j][k];
        atmp[i+1][j][k]=alp[i][j][k]+alp[i+1][j][k];
        fn  [i  ][j][k]+=2;
        fn  [i+1][j][k]+=2;
      }

  /* j-direction */
  for(int j=sj()-1; j<=ej(); j++)
    for_ik(i,k)
      if((clr[i][j][k]-phisurf)*(clr[i][j+1][k]-phisurf)<=0.0){
        atmp[i][j  ][k]=alp[i][j][k]+alp[i][j+1][k];
        atmp[i][j+1][k]=alp[i][j][k]+alp[i][j+1][k];
        fn  [i][j  ][k]+=2;
        fn  [i][j+1][k]+=2;
      }

  /* k-direction */
  for(int k=sk()-1; k<=ek(); k++)
    for_ij(i,j)
      if((clr[i][j][k]-phisurf)*(clr[i][j][k+1]-phisurf)<=0.0){
        atmp[i][j][k  ]=alp[i][j][k]+alp[i][j][k+1];
        atmp[i][j][k+1]=alp[i][j][k]+alp[i][j][k+1];
        fn  [i][j][k  ]+=2;
        fn  [i][j][k+1]+=2;
      }

  for_ijk(i,j,k)
    if(approx(stmp[i][j][k],1.0)){
      alp[i][j][k]=atmp[i][j][k]/fn[i][j][k];
      aindex_max=std::max(alp[i][j][k],aindex_max);
    }

  insert_bc_alp(alp);
  ib_ext_scalar(alp);
  alp.exchange();

  boil::cart.max_real(&aindex_max);
  boil::oout<<"cipcsl2:redist,aindex_max= "<<aindex_max<<"\n";

  if(aindex_max==0.0) return;

  /*-------------+
  |  advect alp  |
  +-------------*/
  for_aijk(i,j,k){
    //if(approx(fabs(iflag[i][j][k]),nlayer)){
    if(abs(iflag[i][j][k])==nlayer){
      stmp[i][j][k]=2.0;
    //} else if(int(fabs(iflag[i][j][k]))>nlayer){
    } else if(abs(iflag[i][j][k])>nlayer){
      stmp[i][j][k]=3.0;
    }
  }

  const real dtau = dxmin;
  real flux,isgn;
  int mmax=4;
  //int mmax=8;
  if(ialpcal==0)
    mmax=20;

  for(int mstep=1; mstep<=mmax; mstep++){
    /* right hand side */
    for_ijk(i,j,k) {
      if(stmp[i][j][k]==0.0){
        isgn=copysign(1.0,dist[i][j][k]);
        flux=-isgn*nx[i][j][k]*dSx(i,j,k);
        atmp[i][j][k] =-(0.5*(alp[i][j][k]+alp[i-1][j][k])*flux
                        +0.5*(alp[i][j][k]-alp[i-1][j][k])*fabs(flux));
        flux= isgn*nx[i][j][k]*dSx(i,j,k);
        atmp[i][j][k]-=0.5*(alp[i][j][k]+alp[i+1][j][k])*flux
                      +0.5*(alp[i][j][k]-alp[i+1][j][k])*fabs(flux);
        flux=-isgn*ny[i][j][k]*dSy(i,j,k);
        atmp[i][j][k]-=0.5*(alp[i][j][k]+alp[i][j-1][k])*flux
                      +0.5*(alp[i][j][k]-alp[i][j-1][k])*fabs(flux);
        flux= isgn*ny[i][j][k]*dSy(i,j,k);
        atmp[i][j][k]-=0.5*(alp[i][j][k]+alp[i][j+1][k])*flux
                      +0.5*(alp[i][j][k]-alp[i][j+1][k])*fabs(flux);
        flux= -isgn*nz[i][j][k]*dSz(i,j,k);
        atmp[i][j][k]-=0.5*(alp[i][j][k]+alp[i][j][k-1])*flux
                      +0.5*(alp[i][j][k]-alp[i][j][k-1])*fabs(flux);
        flux= isgn*nz[i][j][k]*dSz(i,j,k);
        atmp[i][j][k]-=0.5*(alp[i][j][k]+alp[i][j][k+1])*flux
                      +0.5*(alp[i][j][k]-alp[i][j][k+1])*fabs(flux);
      }
    }
    for_aijk(i,j,k)
      fn[i][j][k]=0.0;

    for(int it=1; it<=4; it++){
      int ist,ied,iinc;
      int jst,jed,jinc;
      int kst,ked,kinc;
      if(it%2==0){ist=si();ied=ei();iinc=1;}else{ist=ei();ied=si();iinc=-1;}
      if(it%2==0){jst=sj();jed=ej();jinc=1;}else{jst=ej();jed=sj();jinc=-1;}
      if(it%2==0){kst=sk();ked=ek();kinc=1;}else{kst=ek();ked=sk();kinc=-1;}
      for(int i=ist; i<=ied; i+=iinc){
      for(int j=jst; j<=jed; j+=jinc){
      for(int k=kst; k<=ked; k+=kinc){
        if(stmp[i][j][k]==0.0){
          isgn=copysign(1.0,dist[i][j][k]);
          real diag,drhs;
          diag =dV(i,j,k)/dtau;
          drhs = 0.0;

          flux= -isgn*nx[i][j][k]*dSx(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs =-0.5*(flux-fabs(flux))*fn[i-1][j][k];

          flux= isgn*nx[i][j][k]*dSx(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs-=0.5*(flux-fabs(flux))*fn[i+1][j][k];

          flux=-isgn*ny[i][j][k]*dSy(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs-=0.5*(flux-fabs(flux))*fn[i][j-1][k];

          flux= isgn*ny[i][j][k]*dSy(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs-=0.5*(flux-fabs(flux))*fn[i][j+1][k];

          flux=-isgn*nz[i][j][k]*dSz(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs-=0.5*(flux-fabs(flux))*fn[i][j][k-1];

          flux= isgn*nz[i][j][k]*dSz(i,j,k);
          diag+=0.5*(flux+fabs(flux));
          drhs-=0.5*(flux-fabs(flux))*fn[i][j][k+1];

          fn[i][j][k] = (atmp[i][j][k]+drhs)/diag;
        } else if(approx(stmp[i][j][k],2.0)) {
          real dalpsum=0.0;
          int isum=0;
          if(stmp[i-1][j][k]==0.0){
            dalpsum+=fn[i-1][j][k];
            isum++;
          }
          if(stmp[i+1][j][k]==0.0){
            dalpsum+=fn[i+1][j][k];
            isum++;
          }
          if(stmp[i][j-1][k]==0.0){
            dalpsum+=fn[i][j-1][k];
            isum++;
          }
          if(stmp[i][j+1][k]==0.0){
            dalpsum+=fn[i][j+1][k];
            isum++;
          }
          if(stmp[i][j][k-1]==0.0){
            dalpsum+=fn[i][j][k-1];
            isum++;
          }
          if(stmp[i][j][k+1]==0.0){
            dalpsum+=fn[i][j][k+1];
            isum++;
          }
          fn[i][j][k]=dalpsum/real(isum+1.0e-12);
        }
      }
    }}}

    /* update alp */
    for_ijk(i,j,k) {
      if((stmp[i][j][k]==0.0) || approx(stmp[i][j][k],2.0)){
        alp[i][j][k] += fn[i][j][k];
      } else if(approx(stmp[i][j][k],3.0)) {
        alp[i][j][k]=0.0;
      }
    }
    insert_bc_alp(alp);
    ib_ext_scalar(alp);
    alp.exchange();
  }
  ialpcal=1;

  /*---------------------------+
  | special treatment for wall |
  +---------------------------*/
#ifdef WALL
  /* flag of boundary conditions for stmp & wall */
  bool iminw, imaxw, jminw, jmaxw, kminw, kmaxw;
  iminw = imaxw = jminw = jmaxw = kminw = kmaxw = false;
  for( int b=0; b<phi.bc().count(); b++ ) {
    if( phi.bc().type_decomp(b) ) continue;
    if( phi.bc().type(b)==BndType::wall() ){
      if(phi.bc().direction(b)==Dir::imin()) iminw=true;
      if(phi.bc().direction(b)==Dir::imax()) imaxw=true;
      if(phi.bc().direction(b)==Dir::jmin()) jminw=true;
      if(phi.bc().direction(b)==Dir::jmax()) jmaxw=true;
      if(phi.bc().direction(b)==Dir::kmin()) kminw=true;
      if(phi.bc().direction(b)==Dir::kmax()) kmaxw=true;
    }
  }
  if(iminw){
    int i=si();
    for_jk(j,k) {
      alp[i][j][k]=2.0;
    }
  }
  if(imaxw){
    int i=ei();
    for_jk(j,k) {
      alp[i][j][k]=2.0;
    }
  }
  if(jminw){
    int j=sj();
    for_ik(i,k) {
      alp[i][j][k]=2.0;
    }
  }
  if(jmaxw){
    int j=ej();
    for_ik(i,k) {
      alp[i][j][k]=2.0;
    }
  }
  if(kminw){
    int k=sk();
    for_ij(i,j) {
      alp[i][j][k]=2.0;
    }
  }
  if(kmaxw){
    int k=ek();
    for_ij(i,j) {
      alp[i][j][k]=2.0;
    }
  }
#endif

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    alp[i][j][k]=2.0;
  }
  ib_ext_scalar(alp);
#endif
  alp.exchange();

#ifdef DEBUG
  boil::oout<<"cipcsl2_set_alp:end\n";
#endif
}

  return;
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_set_alp.cpp,v 1.8 2015/05/05 15:12:01 sato Exp $'/
+-----------------------------------------------------------------------------*/
