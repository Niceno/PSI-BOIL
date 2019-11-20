#include "enthalpyfd.h"

using namespace std;

/***************************************************************************//**
*  \brief Adds diffusion to right hand side.
*******************************************************************************/
void EnthalpyFD::update_ftif(const real TS0, const real TSm, const bool nst,
                              const Scalar * diff_eddy) {

  tifmodel.tint_field(nst); 

  /* get time stepping coefficient */
  real tscn = diff_ts.N();
  assert( tscn > 0.0 );

  if( !solid() ) {
    /*------------------------+ 
    |  x direction (w and e)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc,cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrw = std::min(std::max((*clr)[i-1][j][k],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i  ][j][k],0.0),1.0);
      const real clre = std::min(std::max((*clr)[i+1][j][k],0.0),1.0);
      const real pc = phi[i][j][k];
      real xm,xp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      //if((clrw-clrsurf)*(clrc-clrsurf)>=0){
      if(!Interface(-1,Comp::i(),i,j,k)){
        xm=phi.dxw(i);
        pm=phi[i-1][j][k];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clrw-clrc),epsl);
          xm=frac*phi.dxw(i);
          pm=Tint(-1,Comp::i(),frac,i,j,k);
        } else {
          xm = std::max(epsl*phi.dxw(i),distance_x(i,j,k,-1,pm));
        }
        aflagm=1.0;
      }
      //if((clrc-clrsurf)*(clre-clrsurf)>=0){
      if(!Interface(+1,Comp::i(),i,j,k)){
        xp=phi.dxe(i);
        pp=phi[i+1][j][k];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clre-clrc),epsl);
          xp=frac*phi.dxe(i);
          pp=Tint(+1,Comp::i(),frac,i,j,k);
        } else {
          xp = std::max(epsl*phi.dxe(i),distance_x(i,j,k,+1,pp));
        }
        aflagp=1.0;
      }
      real cxm = coef_x_m(xm,xp,phi.xc(i));
      real cxp = coef_x_p(xm,xp,phi.xc(i));

      /* add implicit part of the diffusion term */
      const real Awi = tscn * lc * vol * cxm * aflagm;
      const real Aei = tscn * lc * vol * cxp * aflagp;
      ftif[i][j][k] = Awi*pm + Aei*pp;// - ftifold[i][j][k];
    }

    /*------------------------+ 
    |  y direction (s and n)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc,cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrs = std::min(std::max((*clr)[i][j-1][k],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j  ][k],0.0),1.0);
      const real clrn = std::min(std::max((*clr)[i][j+1][k],0.0),1.0);
      const real pc = phi[i][j][k];
      real ym,yp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      //if((clrs-clrsurf)*(clrc-clrsurf)>=0){
      if(!Interface(-1,Comp::j(),i,j,k)){
        ym=phi.dys(j);
        pm=phi[i][j-1][k];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clrs-clrc),epsl);
          ym=frac*phi.dys(j);
          pm=Tint(-1,Comp::j(),frac,i,j,k);
        } else {
          ym = std::max(epsl*phi.dys(j),distance_y(i,j,k,-1,pm));
        }
        aflagm=1.0;
      }
      //if((clrc-clrsurf)*(clrn-clrsurf)>=0){
      if(!Interface(+1,Comp::j(),i,j,k)){
        yp=phi.dyn(j);
        pp=phi[i][j+1][k];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clrn-clrc),epsl);
          yp=frac*phi.dyn(j);
          pp=Tint(+1,Comp::j(),frac,i,j,k);
        } else {
          yp = std::max(epsl*phi.dyn(j),distance_y(i,j,k,+1,pp));
        }
        aflagp=1.0;
      }
      real cym = coef_y_m(ym,yp,phi.yc(j));
      real cyp = coef_y_p(ym,yp,phi.yc(j));

      const real Asi = tscn * lc * vol * cym * aflagm;
      const real Ani = tscn * lc * vol * cyp * aflagp;
      ftif[i][j][k] += Asi*pm + Ani*pp;
    }

    /*------------------------+ 
    |  z direction (b and t)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc,cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrb = std::min(std::max((*clr)[i][j][k-1],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j][k  ],0.0),1.0);
      const real clrt = std::min(std::max((*clr)[i][j][k+1],0.0),1.0);
      const real pc = phi[i][j][k];
      real zm,zp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      //if((clrb-clrsurf)*(clrc-clrsurf)>=0){
      if(!Interface(-1,Comp::k(),i,j,k)){
        zm=phi.dzb(k);
        pm=phi[i][j][k-1];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clrb-clrc),epsl);
          zm=frac*phi.dzb(k);
          pm=Tint(-1,Comp::k(),frac,i,j,k);
        } else {
          zm = std::max(epsl*phi.dzb(k),distance_z(i,j,k,-1,pm));
        }
        aflagm=1.0;
      }
      //if((clrc-clrsurf)*(clrt-clrsurf)>=0){
      if(!Interface(+1,Comp::k(),i,j,k)){
        zp=phi.dzt(k);
        pp=phi[i][j][k+1];
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clrc)/(clrt-clrc),epsl);
          zp=frac*phi.dzt(k);
          pp=Tint(+1,Comp::k(),frac,i,j,k);
        } else {
          zp = std::max(epsl*phi.dzt(k),distance_z(i,j,k,+1,pp));
        }
        aflagp=1.0;
      }
      real czm = coef_z_m(zm,zp,phi.zc(k));
      real czp = coef_z_p(zm,zp,phi.zc(k));

      const real Abi = tscn * lc * vol * czm * aflagm;
      const real Ati = tscn * lc * vol * czp * aflagp;
      ftif[i][j][k] += Abi*pm + Ati*pp;
    }
    // need to add here immersed boundary without solid !!!

  /*------------------------------------------+
  |  features conduction through solid parts  |
  +------------------------------------------*/
  } else {

    for_m(m){
      int ii,jj,kk;
      ii=jj=kk=0;
      if(m==Comp::i()){
        ii=1;
      } else if(m==Comp::j()){
        jj=1;
      } else {
        kk=1;
      }

      for_ijk(i,j,k) {

        const real vol = phi.dV(i,j,k);
        real lc, xm, xp;

        bool onm, onc, onp, ofm, ofc, ofp; // on & off
        real lsm, lsc, lsp; // lambda
        real clm, clc, clp; // color function
        real dxm, dxp, fdm, fdp, fdms, fdps;
        real edm, edc, edp; // eddy viscosity
        real pm, pc, pp;
        real am, ac, ap;
        real tm, tc, tp;
        real aflagm, aflagp;
        real aream, areap;
        real pos0;
        coef_gen coef_m, coef_p;

        onm=dom->ibody().on (i-ii,j-jj,k-kk);
        onc=dom->ibody().on (i   ,j   ,k   );
        onp=dom->ibody().on (i+ii,j+jj,k+kk);
        ofm=dom->ibody().off(i-ii,j-jj,k-kk);
        ofc=dom->ibody().off(i   ,j   ,k   );
        ofp=dom->ibody().off(i+ii,j+jj,k+kk);
        lsm=solid()->lambda (i-ii,j-jj,k-kk);
        lsc=solid()->lambda (i   ,j   ,k   );
        lsp=solid()->lambda (i+ii,j+jj,k+kk);
        clm=(*clr)[i-ii][j-jj][k-kk];
        clc=(*clr)[i   ][j   ][k   ];
        clp=(*clr)[i+ii][j+jj][k+kk];
	if(m==Comp::i()){
          dxm=phi.dxw(i);
          dxp=phi.dxe(i);
          pos0=phi.xc(i);
          coef_m = &EnthalpyFD::coef_x_m;
          coef_p = &EnthalpyFD::coef_x_p;
          aream = dSx(Sign::neg(),i,j,k);
          areap = dSx(Sign::pos(),i,j,k);
          fdm=dom->ibody().fdxw(i,j,k);
          fdp=dom->ibody().fdxe(i,j,k);
          fdms=dom->ibody().fdxe(i-1,j,k);
          fdps=dom->ibody().fdxw(i+1,j,k);
        } else if(m==Comp::j()){
          dxm=phi.dys(j);
          dxp=phi.dyn(j);
          pos0=phi.yc(j);
          coef_m = &EnthalpyFD::coef_y_m;
          coef_p = &EnthalpyFD::coef_y_p;
          aream = dSy(Sign::neg(),i,j,k);
          areap = dSy(Sign::pos(),i,j,k);
          fdm=dom->ibody().fdys(i,j,k);
          fdp=dom->ibody().fdyn(i,j,k);
          fdms=dom->ibody().fdyn(i,j-1,k);
          fdps=dom->ibody().fdys(i,j+1,k);
        } else {
          dxm=phi.dzb(k);
          dxp=phi.dzt(k);
          pos0=phi.zc(k);
          coef_m = &EnthalpyFD::coef_z_m;
          coef_p = &EnthalpyFD::coef_z_p;
          aream = dSz(Sign::neg(),i,j,k);
          areap = dSz(Sign::pos(),i,j,k);
          fdm=dom->ibody().fdzb(i,j,k);
          fdp=dom->ibody().fdzt(i,j,k);
          fdms=dom->ibody().fdzt(i,j,k-1);
          fdps=dom->ibody().fdzb(i,j,k+1);
        }
        pm=phi[i-ii][j-jj][k-kk];
        pc=phi[i   ][j   ][k   ];
        pp=phi[i+ii][j+jj][k+kk];
        edm=edc=edp=0.0;
        if(diff_eddy){
          if (onm) edm = (*diff_eddy)[i-ii][j-jj][k-kk];
          if (onc) edc = (*diff_eddy)[i   ][j   ][k   ];
          if (onp) edp = (*diff_eddy)[i+ii][j+jj][k+kk];
        }
        diff_matrix(am, ac, ap
                  , tm, tc, tp
                  , aflagm, aflagp
                  , pos0, coef_m, coef_p
                  , vol, aream, areap
                  , onm, onc, onp, ofm, ofc, ofp
                  , lsm, lsc, lsp
                  , clm, clc, clp
                  , dxm, dxp, fdm, fdp, fdms, fdps
                  , pm, pc, pp
                  , edm, edc, edp
                  , i, j, k, m);
  
        ftif[i][j][k] += tscn * (am*(1.0-aflagm)*tm+ap*(1.0-aflagp)*tp);
      }
    }
  }

  for_ijk(i,j,k) {
    ftif[i][j][k] = TS0*ftif[i][j][k]+TSm*ftifold[i][j][k];
  }

  return;
}