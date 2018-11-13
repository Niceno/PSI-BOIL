#include "enthalpytif.h"
using namespace std;

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
                , const real lfm, const real lfc, const real lfp
                , const real clm, const real clc, const real clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m){
  // i,j,k,m: used for debugging
  real lm, lc, lp;                      // lambda
  aflagm=aflagp=1.0;

  /*--------------------+
  |  material property  |
  +--------------------*/
  // minus
  if(onm){
    if(clm>=0.5){
      lm = lambdal + edm*cpl/rhol/turbP;
    } else {
      lm = lambdav + edm*cpv/rhov/turbP;
    }
  } else {
    lm = lsm;
  }

  // center
  if(onc){
    if(clc>=0.5){
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
    if(clp>=0.5){
      lp = lambdal + edp*cpl/rhol/turbP;
    } else {
      lp = lambdav + edp*cpv/rhov/turbP;
    }
  } else {
    lp = lsp;
  }

  if(ofc){
    // center:solid
    tm = pm; 
    tc = pc; 
    tp = pp; 

    if(ofm && ofp){
      // s-s-s
      /* FDM */
      //am = lc * vol * 2.0 / (dxm*(dxm+dxp));
      //ac = lc * vol * 2.0 / (dxm*dxp);
      //ap = lc * vol * 2.0 / (dxp*(dxm+dxp));
      /* FVM */
      am = 0.5 * (lc + lm) * area / dxm;
      ap = 0.5 * (lc + lp) * area / dxp;
      ac = am + ap;

      //std::cout<<"s-s-s: "<<i<<" "<<j<<" "<<k<<"\n";
    } else if(onm && ofp){
      // f-s-s
      fdm = max(fdm,epsl);
      dxm = dxm * fdm;
      /* FDM */
      //am = lc*vol*2.0/(dxm*(dxm+dxp))*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
      //ac = lc*vol*2.0/(dxm*dxp)
      //   - lc*vol*2.0/(dxm*(dxm+dxp))*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
      //ap = lc*vol*2.0/(dxp*(dxm+dxp));
      /* FVM */
      am = lc * area / dxm * fdm*lm/((1.0-fdm)*lc+fdm*lm);
      ap = 0.5 * (lc + lp) * area / dxp;
      ac = am + ap;
      //std::cout<<"f-s-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
    } else if(ofm && onp){
      // s-s-f
      fdp = max(fdp,epsl);
      dxp = dxp * fdp;
      /* FDM */
      //am = lc*vol*2.0/(dxm*(dxm+dxp));
      //ac = lc*vol*2.0/(dxm*dxp)
      //   - lc*vol*2.0/(dxp*(dxm+dxp))*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
      //ap = lc*vol*2.0/(dxp*(dxm+dxp))*fdp*lp/((1.0-fdp)*lc+fdp*lp);
      /* FVM */
      am = 0.5 * (lc + lm) * area / dxm;
      ap = lc * area / dxp * fdp*lp/((1.0-fdp)*lc+fdp*lp);
      ac = am + ap;
      //std::cout<<"s-s-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
    } else {
      // f-s-f
      std::cout<<"diff_matrix: need to be develop!!!\n";
      std::cout<<"fluid-solid-fluid.\n";
      exit(0);
    }

  } else {
    // center:fluid
    tm = pm; 
    tc = pc; 
    tp = pp; 

    if(onm && onp){
      // f-f-f
      if((clm-0.5)*(clc-0.5)>=0){
        dxm=dxm;
      } else {
        if(!fs) {
          real frac = std::max((0.5-clc)/(clm-clc),epsl);
          dxm=frac*dxm;
          tm = Tint_old(-1,m,frac,i,j,k);
        } else {
          #if 0
          if(m==Comp::i())
            dxm = max(epsl*dxm,distance_x(i,j,k,-1,tm,true,&dxm));
          else if(m==Comp::j())
            dxm = max(epsl*dxm,distance_y(i,j,k,-1,tm,true,&dxm));
          else
            dxm = max(epsl*dxm,distance_z(i,j,k,-1,tm,true,&dxm));
          #else
            boil::oout<<"EnthalpyTif:diff_matrix: Underdevelopment. Exiting."<<boil::endl;
            exit(0);
          #endif
        }
        aflagm=0.0;
      }
      if((clc-0.5)*(clp-0.5)>=0){
        dxp=dxp;
      } else {
        if(!fs) {
          real frac = std::max((0.5-clc)/(clp-clc),epsl);
          dxp=frac*dxp;
          tp = Tint_old(+1,m,frac,i,j,k);
        } else {
          #if 0
          if(m==Comp::i())
            dxp = max(epsl*dxp,distance_x(i,j,k,+1,tp,true,&dxp));
          else if(m==Comp::j())
            dxp = max(epsl*dxp,distance_y(i,j,k,+1,tp,true,&dxp));
          else
            dxp = max(epsl*dxp,distance_z(i,j,k,+1,tp,true,&dxp));
          #else
            boil::oout<<"EnthalpyTif:diff_matrix: Underdevelopment. Exiting."<<boil::endl;
            exit(0);
          #endif
        }
        aflagp=0.0;
#if 0
	std::cout<<"aflagp=0.0: " <<i<<" "<<j<<" "<<k<<"\n";
#endif
      }
      if (aflagm==0.0 || aflagp==0.0) {
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
      //std::cout<<"f-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";

    } else if(ofm && onp){ 

      // s-f-f
      fdm = max(fdm,epsl);
      dxm = dxm * fdm;
      if((clc-0.5)*(clp-0.5)>=0){
        dxp=dxp;
      } else {
        if(!fs) {
          real frac = std::max((0.5-clc)/(clp-clc),epsl);
          dxp=frac*dxp;
          tp = Tint_old(+1,m,frac,i,j,k);
        } else {
          #if 0
          if(m==Comp::i())
            dxp = max(epsl*dxp,distance_x(i,j,k,+1,tp,true,&dxp));
          else if(m==Comp::j())
            dxp = max(epsl*dxp,distance_y(i,j,k,+1,tp,true,&dxp));
          else
            dxp = max(epsl*dxp,distance_z(i,j,k,+1,tp,true,&dxp));
          #else
            boil::oout<<"EnthalpyTif:diff_matrix: Underdevelopment. Exiting."<<boil::endl;
            exit(0);
          #endif
        }
        aflagp=0.0;
      }
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
      //std::cout<<"s-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
 
    } else if(onm && ofp){

      // f-f-s
      fdp = max(fdp,epsl);
      dxp = dxp * fdp;
      if((clm-0.5)*(clc-0.5)>=0){
        dxm=dxm;
      } else {
        if(!fs) {
          real frac = std::max((0.5-clc)/(clm-clc),epsl);
          dxm=frac*dxm;
          tm = Tint_old(-1,m,frac,i,j,k);
        } else {
          #if 0
          if(m==Comp::i())
            dxm = max(epsl*dxm,distance_x(i,j,k,-1,tm,true,&dxm));
          else if(m==Comp::j())
            dxm = max(epsl*dxm,distance_y(i,j,k,-1,tm,true,&dxm));
          else
            dxm = max(epsl*dxm,distance_z(i,j,k,-1,tm,true,&dxm));
          #else
            boil::oout<<"EnthalpyTif:diff_matrix: Underdevelopment. Exiting."<<boil::endl;
            exit(0);
          #endif
        }
        aflagm=0.0;
      }
      if (aflagm==0.0) {
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
      //std::cout<<"f-f-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
 
    } else {

      // s-f-s
      std::cout<<"diff_matrix: need to be develop!!!\n";
      std::cout<<"solid-fluid-solid.\n";
      exit(0);

    }
  }

  return;
}
