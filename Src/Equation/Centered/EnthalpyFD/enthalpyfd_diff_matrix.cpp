#include "enthalpyfd.h"

//#define DEBUG
//#define USE_FDM_SOLID
//#define USE_FDM_FLUID
/* using FDM worsens the performance. Also, for FDM, the correct BndGrid for Ne
   umann boundary condition is extrapolate, rather than wall */

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyFD::diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , real & aflagm, real & aflagp
                , real & sourceterm
                , const real x0, const coef_gen coef_m, const coef_gen coef_p
                , const real vol, const real aream, const real areap
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const int clm, const int clc, const int clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m) {
  real lm, lc, lp; // lambda
  aflagm=aflagp=1.0;

#ifdef DEBUG
  boil::oout<<"start: "<<i<<" "<<j<<" "<<k<<" "<<m<<" | "<<onm<<" "<<onc<<" "<<onp<<boil::endl;
#endif

  /*--------------------+
  |  material property  |
  +--------------------*/
  // minus
  if(onm){
    if(clm>0){
      lm = lambdal + edm*cpl/rhol/turbP;
    } else {
      lm = lambdav + edm*cpv/rhov/turbP;
    }
  } else {
    lm = lsm;
  }

  // center
  if(onc){
    if(clc>0){
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
    if(clp>0){
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
#ifdef USE_FDM_SOLID
      /* FDM */
      am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
      ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
      ac = am + ap;
#else
      /* FVM */
      am = 0.5 * (lc + lm) * aream / dxm;
      ap = 0.5 * (lc + lp) * areap / dxp;
      ac = am + ap;
#endif
#ifdef DEBUG
      std::cout<<"s-s-s: "<<i<<" "<<j<<" "<<k<<"\n";
#endif
    } else if(onm && ofp){

      real fact;
      real resistm;

      /* f-s-s */
      if(interface(Sign::neg(),m,i,j,k)) {
        /* removed from system matrix */
        aflagm = 0.0;
        /* dxm,fdm are corrected to account for interface position */
        /* tm is changed to the saturation temperature */
        fdm *= dxm;
        dxm = distance_int(Sign::neg(),m,i,j,k,tm);
        real dxfull = dxm;
        fdm /= dxm; 
        fdm = std::max(fdm,epsl);
        /* dxm is corrected */
        dxm = dxm * fdm;
        /* lambdaf is inverted */
        if(clm>0) {
          lm = lambdav + edm*cpv/rhov/turbP;
        } else {
          lm = lambdal + edm*cpl/rhol/turbP;
        }
        fact = resistance_multiplier(dxm,dxfull-dxm,lc,lm,
                                     htwallmodel->near_wall_resist);
        resistm = (dxfull-dxm)/lm + htwallmodel->near_wall_resist;
      } else {
        fdm = std::max(fdm,epsl);
        resistm = (1.0-fdm)*dxm/lm;
        dxm = dxm * fdm;
        fact = fdm*lm/(fdm*lm+(1.0-fdm)*lc);
      }

#ifdef USE_FDM_SOLID
      /* FDM */
      am = lc*vol*(this->*coef_m)(dxm,dxp,x0)*fact;
      ap = lc*vol*(this->*coef_p)(dxm,dxp,x0);
      ac = am + ap;
#else
      /* FVM */
      am = lc * aream / dxm * fact;
      ap = 0.5 * (lc + lp) * areap / dxp;
      ac = am + ap;
#endif
      /* dirac source term */
      sourceterm = htwallmodel->dirac_wall_source * am * resistm;

#ifdef DEBUG
      std::cout<<"f-s-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(ofm && onp){

      real fact;
      real resistp;

      /* s-s-f */
      if(interface(Sign::pos(),m,i,j,k)) {
        /* removed from system matrix */
        aflagp = 0.0;
        /* dxp,fdp are corrected to account for interface position */
        /* tp is changed to the saturation temperature */
        fdp *= dxp;
        dxp = distance_int(Sign::pos(),m,i,j,k,tp);
        real dxfull = dxp;
        fdp /= dxp; 
        fdp = std::max(fdp,epsl);
        /* dxp is corrected */
        dxp = dxp * fdp;
        /* lambdaf is inverted */
        if(clp>0){
          lp = lambdav + edp*cpv/rhov/turbP;
        } else {
          lp = lambdal + edp*cpl/rhol/turbP;
        }
        fact = resistance_multiplier(dxp,dxfull-dxp,lc,lp,
                                     htwallmodel->near_wall_resist);
        resistp = (dxfull-dxp)/lp + htwallmodel->near_wall_resist;
#ifdef DEBUG
        boil::oout<<"fact: "<<i<<" "<<j<<" "<<k<<" "
                  << fact<<" "<<fdp*lp/((1.0-fdp)*lc+fdp*lp)<<boil::endl;
#endif
      } else {
        fdp = std::max(fdp,epsl);
        resistp = (1.0-fdp)*dxp/lp;
        dxp = dxp * fdp;
        fact = fdp*lp/((1.0-fdp)*lc+fdp*lp);
      }

#ifdef USE_FDM_SOLID
      /* FDM */
      am = lc*vol*(this->*coef_m)(dxm,dxp,x0);
      ap = lc*vol*(this->*coef_p)(dxm,dxp,x0)*fact;
      ac = am + ap;
#else
      /* FVM */
      am = 0.5 * (lc + lm) * aream / dxm;
      ap = lc * areap / dxp * fact;
      ac = am + ap;
#endif
      /* dirac source term */
      sourceterm = htwallmodel->dirac_wall_source * ap * resistp;

#ifdef DEBUG
      std::cout<<"s-s-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else {
      /* f-s-f */
      std::cout<<"EnthFD::diff_matrix: Underdevelopment!\n";
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
      if(!interface(Sign::neg(),m,i,j,k)){
        dxm=dxm;
      } else {
        dxm = distance_int(Sign::neg(),m,i,j,k,tm);
        aflagm=0.0;
#ifdef DEBUG
	std::cout<<"aflagm=0.0: " <<i<<" "<<j<<" "<<k<<" "<<clc<<"\n";
#endif
      }
      if(!interface(Sign::pos(),m,i,j,k)){
        dxp=dxp;
      } else {
        dxp = distance_int(Sign::pos(),m,i,j,k,tp);
        aflagp=0.0;
#ifdef DEBUG
	std::cout<<"aflagp=0.0: " <<i<<" "<<j<<" "<<k<<" "<<clc<<"\n";
#endif
      }
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
      /* FDM */
      am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
      ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
      ac = am + ap;
#else
      //if(aflagm==0.0 || aflagp==0.0) {
      if(fabs(clc)<3) {
        /* FDM */
        am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
        ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
      } else { 
        /* FVM */
        am = 0.5 * (lc + lm) * aream / dxm;
        ap = 0.5 * (lc + lp) * areap / dxp;
        ac = am + ap;
      }
#endif
#ifdef DEBUG
      std::cout<<"f-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(ofm && onp){ 
      /* s-f-f */
      if(!interface(Sign::pos(),m,i,j,k)){
        dxp=dxp;
      } else {
        dxp = distance_int(Sign::pos(),m,i,j,k,tp);
        aflagp=0.0;
      }
      if(interface(Sign::neg(),m,i,j,k)) {
        /* removed from system matrix */
        aflagm = 0.0;
        /* dxm is corrected to account for interface position */
        /* tm is changed to the saturation temperature */
        dxm = distance_int(Sign::neg(),m,i,j,k,tm);
        /* FDM */
        am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
        ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
      } else {
        fdm = std::max(fdm,epsl);
        real resistm = (1.0-fdm)*dxm/lm;
        dxm = dxm * fdm;
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*(this->*coef_m)(dxm,dxp,x0)*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
        ac = lc*vol*((this->*coef_m)(dxm,dxp,x0)+(this->*coef_p)(dxm,dxp,x0))
           - lc*vol*(this->*coef_m)(dxm,dxp,x0)*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
        ap = lc*vol*(this->*coef_p)(dxm,dxp,x0);
#else
        //if(aflagp==0.0) {
        if(fabs(clc)<3) {
          /* FDM */
          am = lc*vol*(this->*coef_m)(dxm,dxp,x0)*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
          ac = lc*vol*((this->*coef_m)(dxm,dxp,x0)+(this->*coef_p)(dxm,dxp,x0))
             - lc*vol*(this->*coef_m)(dxm,dxp,x0)*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
          ap = lc*vol*(this->*coef_p)(dxm,dxp,x0);
        } else {
          /* FVM */
          am = lc * aream / dxm * fdm * lm / (fdm*lm+(1.0-fdm)*lc);
          ap = 0.5 * (lc + lp) * areap / dxp;
          ac = am + ap;
        }
#endif
        /* dirac source term */
        sourceterm = htwallmodel->dirac_wall_source * am * resistm;
      }
#ifdef DEBUG
      std::cout<<"s-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(onm && ofp){
      /* f-f-s */
      if(!interface(Sign::neg(),m,i,j,k)){
        dxm=dxm;
      } else {
        dxm = distance_int(Sign::neg(),m,i,j,k,tm);
        aflagm=0.0;
      }
      if(interface(Sign::pos(),m,i,j,k)) {
        /* removed from system matrix */
        aflagp = 0.0;
        /* dxp is corrected to account for interface position */
        /* tp is changed to the saturation temperature */
        dxp = distance_int(Sign::pos(),m,i,j,k,tp);
        /* FDM */
        am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
        ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
      } else {
        fdp = std::max(fdp,epsl);
        real resistp = (1.0-fdp)*dxp/lp;
        dxp = dxp * fdp;
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*(this->*coef_m)(dxm,dxp,x0);
        ac = lc*vol*((this->*coef_m)(dxm,dxp,x0)+(this->*coef_p)(dxm,dxp,x0))
           - lc*vol*(this->*coef_p)(dxm,dxp,x0)*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
        ap = lc*vol*(this->*coef_p)(dxm,dxp,x0)*fdp*lp/((1.0-fdp)*lc+fdp*lp);
#else
        //if(aflagm==0.0) {
        if(fabs(clc)<3) {
          /* FVM */
          am = lc*vol*(this->*coef_m)(dxm,dxp,x0);
          ac = lc*vol*((this->*coef_m)(dxm,dxp,x0)+(this->*coef_p)(dxm,dxp,x0))
             - lc*vol*(this->*coef_p)(dxm,dxp,x0)*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
          ap = lc*vol*(this->*coef_p)(dxm,dxp,x0)*fdp*lp/((1.0-fdp)*lc+fdp*lp);
        } else {
          /* FVM */
          am = 0.5 * (lc + lm) * aream / dxm;
          ap = lc * areap / dxp * fdp * lp / (fdp*lp+(1.0-fdp)*lc);
          ac = am + ap;
        }
#endif
        /* dirac source term */
        sourceterm = htwallmodel->dirac_wall_source * ap * resistp;
      }
#ifdef DEBUG
      std::cout<<"f-f-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif 
    } else {
      /* s-f-s */
      std::cout<<"EnthFD::diff_matrix: Underdevelopment!\n";
      std::cout<<"solid-fluid-solid.\n";
      exit(0);
    }
  }

#ifdef DEBUG
  boil::oout<<"end: "<<i<<" "<<j<<" "<<k<<" "<<m<<boil::endl;
#endif

  return;
}

/***************************************************************************//**
* Effect of flux continuity 
*******************************************************************************/
inline real EnthalpyFD::resistance_multiplier(const real dx1, const real dx2,
                                              const real l1, const real l2,
                                              const real resistplus) const {
  return dx1/l1 / (dx1/l1 + dx2/l2 + resistplus);
}
