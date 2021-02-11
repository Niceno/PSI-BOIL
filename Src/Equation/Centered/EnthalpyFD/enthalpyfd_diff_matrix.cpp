#include "enthalpyfd.h"

#include "def.h"
//#define DEBUG

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyFD::diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , bool & aflagm, bool & aflagp
                , real & sourceterm
                , const real x0, const coef_gen coef_m, const coef_gen coef_p
                , const real vol, const real aream, const real areap
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const real lvm, const real lvc, const real lvp
                , const real llm, const real llc, const real llp
                , const int clm, const int clc, const int clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , const int i, const int j, const int k, const Comp m) {
  real lm, lc, lp; // lambda
  aflagm=aflagp=true;

#ifdef DEBUG
  boil::oout<<"start: "<<i<<" "<<j<<" "<<k<<" "<<m<<" | "<<onm<<" "<<onc<<" "<<onp<<boil::endl;
#endif

  /*--------------------+
  |  material property  |
  +--------------------*/
  // minus
  if(onm){
    if(clm>0){
      lm = llm;
    } else {
      lm = lvm;
    }
  } else {
    lm = lsm;
  }

  // center
  if(onc){
    if(clc>0){
      lc = llc;
    } else {
      lc = lvc;
    }
  } else {
    lc = lsc;
    fdm = 1.0-fdms;
    fdp = 1.0-fdps;
  }

  // plus
  if(onp){
    if(clp>0){
      lp = llp;
    } else {
      lp = lvp;
    }
  } else {
    lp = lsp;
  }

  /*------------------+
  |  center in solid  |
  +------------------*/
  if(ofc){
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

      /* f-s-s */
      real fact;
      real resist_f, resist_s;
      Sign dummy; /* dummy sign */

      /* fluid resistance */
      if(cht.interface(Sign::neg(),m,i,j,k)) {
        /* removed from system matrix */
        aflagm = false;
        /* tm is changed to the saturation temperature */
        real dxfull = cht.distance_int(Sign::neg(),m,i,j,k,tm,dummy,
                                       ResistEval::no,Old::no);
        dxm = dxm * fdm;
        /* lambdaf is inverted */
        if(clm>0) {
          lm = lvm;
        } else {
          lm = llm;
        }
        resist_f = (dxfull-dxm)/lm;
        resist_s = dxm/lc;
        resist_f += cht.wall_resistance(i-(m==Comp::i()),
                                        j-(m==Comp::j()),
                                        k-(m==Comp::k()));
      } else {
        resist_f = (1.0-fdm)*dxm/lm;
        dxm = dxm * fdm;
        resist_s = dxm/lc;
        resist_f += cht.wall_resistance(i-(m==Comp::i()),
                                        j-(m==Comp::j()),
                                        k-(m==Comp::k()));
      }
      fact = resistance_multiplier(resist_s,resist_f);

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
      sourceterm = cht.dirac_wall_source(i-(m==Comp::i()),
                                         j-(m==Comp::j()),
                                         k-(m==Comp::k())) * am * resist_f;

#ifdef DEBUG
      std::cout<<"f-s-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(ofm && onp){

      /* s-s-f */
      real fact;
      real resist_f, resist_s;
      Sign dummy; /* dummy sign */

      /* fluid resistance */
      if(cht.interface(Sign::pos(),m,i,j,k)) {
        /* removed from system matrix */
        aflagp = false;
        /* tp is changed to the saturation temperature */
        real dxfull = cht.distance_int(Sign::pos(),m,i,j,k,tp,dummy,
                                       ResistEval::no,Old::no);
        dxp = dxp * fdp;
        /* lambdaf is inverted */
        if(clp>0){
          lp = lvp;
        } else {
          lp = llp;
        }
        resist_f = (dxfull-dxp)/lp;
        resist_s = dxp/lc;
        resist_f += cht.wall_resistance(i+(m==Comp::i()),
                                        j+(m==Comp::j()),
                                        k+(m==Comp::k()));
      } else {
        resist_f = (1.0-fdp)*dxp/lp;
        dxp = dxp * fdp;
        resist_s = dxp/lc;
        resist_f += cht.wall_resistance(i+(m==Comp::i()),
                                        j+(m==Comp::j()),
                                        k+(m==Comp::k()));
      }
      fact = resistance_multiplier(resist_s,resist_f);
#ifdef DEBUG
      boil::oout<<"fact: "<<i<<" "<<j<<" "<<k<<" "
                << fact<<" "<<fdp*lp/((1.0-fdp)*lc+fdp*lp)<<boil::endl;
#endif

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
      sourceterm = cht.dirac_wall_source(i+(m==Comp::i()),
                                         j+(m==Comp::j()),
                                         k+(m==Comp::k())) * ap * resist_f;

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
    Sign cell_marker;

    if(onm && onp){
      /* f-f-f */
      if(!cht.interface(Sign::neg(),m,i,j,k)){
        dxm=dxm;
      } else {
        dxm = cht.distance_int(Sign::neg(),m,i,j,k,tm,cell_marker,
                               ResistEval::no,Old::no);
        aflagm=false;
#ifdef DEBUG
	std::cout<<"aflagm=false: " <<i<<" "<<j<<" "<<k<<" "<<clc<<"\n";
#endif
      }
      if(!cht.interface(Sign::pos(),m,i,j,k)){
        dxp=dxp;
      } else {
        dxp = cht.distance_int(Sign::pos(),m,i,j,k,tp,cell_marker,
                               ResistEval::no,Old::no);
        aflagp=false;
#ifdef DEBUG
	std::cout<<"aflagp=false: " <<i<<" "<<j<<" "<<k<<" "<<clc<<"\n";
#endif
      }
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
      /* FDM */
      am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
      ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
      ac = am + ap;
#else
      //if(!aflagm || !aflagp) {
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
      if(!cht.interface(Sign::pos(),m,i,j,k)){
        dxp=dxp;
      } else {
        dxp = cht.distance_int(Sign::pos(),m,i,j,k,tp,cell_marker,
                               ResistEval::no,Old::no);
        aflagp=false;
      }
      if(cht.interface(Sign::neg(),m,i,j,k)) {
        /* removed from system matrix */
        aflagm = false;
        /* dxm is corrected to account for interface position */
        /* tm is changed to the saturation temperature */
        dxm = cht.distance_int(Sign::neg(),m,i,j,k,tm,cell_marker,
                               ResistEval::no,Old::no);
        /* FDM */
        am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
        ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
      } else {
        real resist_s = (1.0-fdm)*dxm/lm + cht.wall_resistance(i,j,k);
        dxm = dxm * fdm;
        real resist_f = dxm/lc;
        real fact = resistance_multiplier(resist_f,resist_s);
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*(this->*coef_m)(dxm,dxp,x0)*fact;
        ap = lc*vol*(this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
#else
        //if(!aflagp) {
        if(fabs(clc)<3) {
          /* FDM */
          am = lc*vol*(this->*coef_m)(dxm,dxp,x0)*fact;
          ap = lc*vol*(this->*coef_p)(dxm,dxp,x0);
          ac = am + ap;
        } else {
          /* FVM */
          am = lc * aream / dxm * fact;
          ap = 0.5 * (lc + lp) * areap / dxp;
          ac = am + ap;
        }
#endif
        /* dirac source term */
        sourceterm = cht.dirac_wall_source(i,j,k) * am * resist_s;
      }
#ifdef DEBUG
      std::cout<<"s-f-f: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
#endif
    } else if(onm && ofp){
      /* f-f-s */
      if(!cht.interface(Sign::neg(),m,i,j,k)){
        dxm=dxm;
      } else {
        dxm = cht.distance_int(Sign::neg(),m,i,j,k,tm,cell_marker,
                               ResistEval::no,Old::no);
        aflagm=false;
      }
      if(cht.interface(Sign::pos(),m,i,j,k)) {
        /* removed from system matrix */
        aflagp = false;
        /* dxp is corrected to account for interface position */
        /* tp is changed to the saturation temperature */
        dxp = cht.distance_int(Sign::pos(),m,i,j,k,tp,cell_marker,
                               ResistEval::no,Old::no);
        /* FDM */
        am = lc * vol * (this->*coef_m)(dxm,dxp,x0);
        ap = lc * vol * (this->*coef_p)(dxm,dxp,x0);
        ac = am + ap;
      } else {
        real resist_s = (1.0-fdp)*dxp/lp + cht.wall_resistance(i,j,k);
        dxp = dxp * fdp;
        real resist_f = dxp/lc;
        real fact = resistance_multiplier(resist_f,resist_s);
#ifdef USE_FDM_FLUID /* now, finite difference is applied always, as in system_diffusive */
        /* FDM */
        am = lc*vol*(this->*coef_m)(dxm,dxp,x0);
        ap = lc*vol*(this->*coef_p)(dxm,dxp,x0)*fact;
        ac = am + ap;
#else
        //if(!aflagm) {
        if(fabs(clc)<3) {
          /* FVM */
          am = lc*vol*(this->*coef_m)(dxm,dxp,x0);
          ap = lc*vol*(this->*coef_p)(dxm,dxp,x0)*fact;
          ac = am + ap;
        } else {
          /* FVM */
          am = 0.5 * (lc + lm) * aream / dxm;
          ap = lc * areap / dxp * fact;
          ac = am + ap;
        }
#endif
        /* dirac source term */
        sourceterm = cht.dirac_wall_source(i,j,k) * ap * resist_s;
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
