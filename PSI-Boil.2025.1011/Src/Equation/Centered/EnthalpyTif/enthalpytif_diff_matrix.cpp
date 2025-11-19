#include "enthalpytif.h"
using namespace std;

#define USE_FDM

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyTIF::diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , real & aflagm, real & aflagp
                , const real vol, const real area
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const real clm, const real clc, const real clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m) {
  // i,j,k,m: used for debugging
  real lm, lc, lp;                      // lambda
  aflagm=aflagp=1.0;

  /*--------------------+
  |  material property  |
  +--------------------*/
  // minus
  if(onm){
    if(clm>=clrsurf){
      lm = lambdal + edm*cpl/rhol/turbP;
    } else {
      lm = lambdav + edm*cpv/rhov/turbP;
    }
  } else {
    lm = lsm;
  }

  // center
  if(onc){
    if(clc>=clrsurf){
      lc = lambdal + edc*cpl/rhol/turbP;
    } else {
      lc = lambdav + edc*cpv/rhov/turbP;
    }
  } else {
    lc = lsc;
    fdm = 1.0-fdms;
    fdp = 1.0-fdps;
  }

  // plus
  if(onp){
    if(clp>=clrsurf){
      lp = lambdal + edp*cpl/rhol/turbP;
    } else {
      lp = lambdav + edp*cpv/rhov/turbP;
    }
  } else {
    lp = lsp;
  }

  /*------------------+
  |  center in solid  |
  +------------------*/
  if(ofc){
    tm = pm; 
    tc = pc; 
    tp = pp; 

    if(ofm && ofp){
      /* s-s-s */
#ifdef USE_FDM
      /* FDM */
      am = lc * vol * 2.0 / (dxm*(dxm+dxp));
      ac = lc * vol * 2.0 / (dxm*dxp);
      ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
#else
      /* FVM */
      am = 0.5 * (lc + lm) * area / dxm;
      ap = 0.5 * (lc + lp) * area / dxp;
      ac = am + ap;
#endif
#if 0
      std::cout<<"s-s-s: "<<i<<" "<<j<<" "<<k<<"\n";
#endif
    } else if(onm && ofp){
      /* f-s-s */
      if(fs&&Interface(-1,m,i,j,k)) {
        /* removed from system matrix */
        aflagm = 0.0;
        /* dxm,fdm are corrected to account for interface position */
        /* tm is changed to the saturation temperature */
        fdm *= dxm;
        if(m==Comp::i()) 
          dxm = std::max(epsl*dxm,distance_x(i,j,k,-1,tm));
        else if(m==Comp::j())
          dxm = std::max(epsl*dxm,distance_y(i,j,k,-1,tm));
        else
          dxm = std::max(epsl*dxm,distance_z(i,j,k,-1,tm));
        fdm /= dxm; 
        fdm = max(fdm,epsl);
        /* dxm is corrected */
        dxm = dxm * fdm;
        /* lambdaf is inverted */
        if(clm>=clrsurf){
          lm = lambdav + edm*cpv/rhov/turbP;
        } else {
          lm = lambdal + edm*cpl/rhol/turbP;
        }
      } else {
        fdm = max(fdm,epsl);
        dxm = dxm * fdm;
      }
#ifdef USE_FDM
      /* FDM */
      am = lc*vol*2.0/(dxm*(dxm+dxp))*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
      ac = lc*vol*2.0/(dxm*dxp)
         - lc*vol*2.0/(dxm*(dxm+dxp))*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
      ap = lc*vol*2.0/(dxp*(dxm+dxp));
#else
      /* FVM */
      am = lc * area / dxm * fdm*lm/((1.0-fdm)*lc+fdm*lm);
      ap = 0.5 * (lc + lp) * area / dxp;
      ac = am + ap;
#endif
#if 0
      std::cout<<"f-s-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(ofm && onp){
      /* s-s-f */
      if(fs&&Interface(+1,m,i,j,k)) {
        /* removed from system matrix */
        aflagp = 0.0;
        /* dxp,fdp are corrected to account for interface position */
        /* tp is changed to the saturation temperature */
        fdp *= dxp;
        if(m==Comp::i()) 
          dxp = std::max(epsl*dxp,distance_x(i,j,k,+1,tp));
        else if(m==Comp::j())
          dxp = std::max(epsl*dxp,distance_y(i,j,k,+1,tp));
        else
          dxp = std::max(epsl*dxp,distance_z(i,j,k,+1,tp));
        fdp /= dxp; 
        fdp = max(fdp,epsl);
        /* dxp is corrected */
        dxp = dxp * fdp;
        /* lambdaf is inverted */
        if(clp>=clrsurf){
          lp = lambdav + edp*cpv/rhov/turbP;
        } else {
          lp = lambdal + edp*cpl/rhol/turbP;
        }
      } else {
        fdp = max(fdp,epsl);
        dxp = dxp * fdp;
      }
#ifdef USE_FDM
      /* FDM */
      am = lc*vol*2.0/(dxm*(dxm+dxp));
      ac = lc*vol*2.0/(dxm*dxp)
         - lc*vol*2.0/(dxp*(dxm+dxp))*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
      ap = lc*vol*2.0/(dxp*(dxm+dxp))*fdp*lp/((1.0-fdp)*lc+fdp*lp);
#else
      /* FVM */
      am = 0.5 * (lc + lm) * area / dxm;
      ap = lc * area / dxp * fdp*lp/((1.0-fdp)*lc+fdp*lp);
      ac = am + ap;
#endif
#if 0
      std::cout<<"s-s-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else {
      /* f-s-f */
      std::cout<<"ETIF::diff_matrix: Underdevelopment!\n";
      std::cout<<"fluid-solid-fluid.\n";
      exit(0);
    }

  /*------------------+
  |  center in fluid  |
  +------------------*/
  } else {
    tm = pm; 
    tc = pc; 
    tp = pp; 

    if(onm && onp){
      /* f-f-f */
      //if((clm-clrsurf)*(clc-clrsurf)>=0){
      if(!Interface(-1,m,i,j,k)){
        dxm=dxm;
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clc)/(clm-clc),epsl);
          dxm=frac*dxm;
          tm = Tint(-1,m,frac,i,j,k);
        } else {
          if(m==Comp::i())
            dxm = std::max(epsl*dxm,distance_x(i,j,k,-1,tm));
          else if(m==Comp::j())
            dxm = std::max(epsl*dxm,distance_y(i,j,k,-1,tm));
          else
            dxm = std::max(epsl*dxm,distance_z(i,j,k,-1,tm));
        }
        aflagm=0.0;
      }
      //if((clc-clrsurf)*(clp-clrsurf)>=0){
      if(!Interface(+1,m,i,j,k)){
        dxp=dxp;
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clc)/(clp-clc),epsl);
          dxp=frac*dxp;
          tp = Tint(+1,m,frac,i,j,k);
        } else {
          if(m==Comp::i())
            dxp = std::max(epsl*dxp,distance_x(i,j,k,+1,tp));
          else if(m==Comp::j())
            dxp = std::max(epsl*dxp,distance_y(i,j,k,+1,tp));
          else
            dxp = std::max(epsl*dxp,distance_z(i,j,k,+1,tp));
        }
        aflagp=0.0;
#if 0
	std::cout<<"aflagp=0.0: " <<i<<" "<<j<<" "<<k<<"\n";
#endif
      }
#ifdef USE_FDM /* now, finite difference is applied always, as in system_diffusive */
      /* FDM */
      am = lc * vol * 2.0 / (dxm*(dxm+dxp));
      ac = lc * vol * 2.0 / (dxm*dxp);
      ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
#else
      if(aflagm==0.0 || aflagp==0.0) {
        /* FDM */
        am = lc * vol * 2.0 / (dxm*(dxm+dxp));
        ac = lc * vol * 2.0 / (dxm*dxp);
        ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
      } else { 
        /* FVM */
        am = 0.5 * (lc + lm) * area / dxm;
        ap = 0.5 * (lc + lp) * area / dxp;
        ac = am + ap;
      }
#endif
#if 0
      std::cout<<"f-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(ofm && onp){ 
      /* s-f-f */
      //if((clc-clrsurf)*(clp-clrsurf)>=0){
      if(!Interface(+1,m,i,j,k)){
        dxp=dxp;
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clc)/(clp-clc),epsl);
          dxp=frac*dxp;
          tp = Tint(+1,m,frac,i,j,k);
        } else {
          if(m==Comp::i())
            dxp = std::max(epsl*dxp,distance_x(i,j,k,+1,tp));
          else if(m==Comp::j())
            dxp = std::max(epsl*dxp,distance_y(i,j,k,+1,tp));
          else
            dxp = std::max(epsl*dxp,distance_z(i,j,k,+1,tp));
        }
        aflagp=0.0;
      }
      if(fs&&Interface(-1,m,i,j,k)) {
        /* removed from system matrix */
        aflagm = 0.0;
        /* dxm is corrected to account for interface position */
        /* tm is changed to the saturation temperature */
        if(m==Comp::i()) 
          dxm = std::max(epsl*dxm,distance_x(i,j,k,-1,tm));
        else if(m==Comp::j())
          dxm = std::max(epsl*dxm,distance_y(i,j,k,-1,tm));
        else
          dxm = std::max(epsl*dxm,distance_z(i,j,k,-1,tm));
        /* FDM */
        am = lc * vol * 2.0 / (dxm*(dxm+dxp));
        ac = lc * vol * 2.0 / (dxm*dxp);
        ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
      } else {
        fdm = max(fdm,epsl);
        dxm = dxm * fdm;
#ifdef USE_FDM /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*2.0/(dxm*(dxm+dxp))*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
        ac = lc*vol*2.0/(dxm*dxp)
           - lc*vol*2.0/(dxm*(dxm+dxp))*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
        ap = lc*vol*2.0/(dxp*(dxm+dxp));
#else
        if(aflagp==0.0) {
          /* FDM */
          am = lc*vol*2.0/(dxm*(dxm+dxp))*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
          ac = lc*vol*2.0/(dxm*dxp)
             - lc*vol*2.0/(dxm*(dxm+dxp))*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
          ap = lc*vol*2.0/(dxp*(dxm+dxp));
        } else {
          /* FVM */
          am = lc * area / dxm * fdm * lm / (fdm*lm+(1.0-fdm)*lc);
          ap = 0.5 * (lc + lp) * area / dxp;
          ac = am + ap;
        }
#endif
      }
#if 0
      std::cout<<"s-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(onm && ofp){
      /* f-f-s */
      //if((clm-clrsurf)*(clc-clrsurf)>=0){
      if(!Interface(-1,m,i,j,k)){
        dxm=dxm;
      } else {
        if(!fs) {
          real frac = std::max((clrsurf-clc)/(clm-clc),epsl);
          dxm=frac*dxm;
          tm = Tint(-1,m,frac,i,j,k);
        } else {
          if(m==Comp::i())
            dxm = std::max(epsl*dxm,distance_x(i,j,k,-1,tm));
          else if(m==Comp::j())
            dxm = std::max(epsl*dxm,distance_y(i,j,k,-1,tm));
          else
            dxm = std::max(epsl*dxm,distance_z(i,j,k,-1,tm));
        }
        aflagm=0.0;
      }
      if(fs&&Interface(+1,m,i,j,k)) {
        /* removed from system matrix */
        aflagp = 0.0;
        /* dxp is corrected to account for interface position */
        /* tp is changed to the saturation temperature */
        if(m==Comp::i())
          dxp = std::max(epsl*dxp,distance_x(i,j,k,+1,tp));
        else if(m==Comp::j())
          dxp = std::max(epsl*dxp,distance_y(i,j,k,+1,tp));
        else
          dxp = std::max(epsl*dxp,distance_z(i,j,k,+1,tp));
        /* FDM */
        am = lc * vol * 2.0 / (dxm*(dxm+dxp));
        ac = lc * vol * 2.0 / (dxm*dxp);
        ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
      } else {
        fdp = max(fdp,epsl);
        dxp = dxp * fdp;
#ifdef USE_FDM /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*2.0/(dxm*(dxm+dxp));
        ac = lc*vol*2.0/(dxm*dxp)
           - lc*vol*2.0/(dxp*(dxm+dxp))*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
        ap = lc*vol*2.0/(dxp*(dxm+dxp))*fdp*lp/((1.0-fdp)*lc+fdp*lp);
#else
        if(aflagm==0.0) {
          /* FVM */
          am = lc*vol*2.0/(dxm*(dxm+dxp));
          ac = lc*vol*2.0/(dxm*dxp)
             - lc*vol*2.0/(dxp*(dxm+dxp))*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
          ap = lc*vol*2.0/(dxp*(dxm+dxp))*fdp*lp/((1.0-fdp)*lc+fdp*lp);
        } else {
          /* FVM */
          am = 0.5 * (lc + lm) * area / dxm;
          ap = lc * area / dxp * fdp * lp / (fdp*lp+(1.0-fdp)*lc);
          ac = am + ap;
        }
#endif
      }
#if 0
      std::cout<<"f-f-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif 
    } else {
      /* s-f-s */
      std::cout<<"ETIF::diff_matrix: Underdevelopment!\n";
      std::cout<<"solid-fluid-solid.\n";
      exit(0);
    }
  }

  return;
}
