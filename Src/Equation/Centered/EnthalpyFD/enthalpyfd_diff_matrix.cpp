#include "enthalpyfd.h"
using namespace boil;

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyFD::diff_matrix(real & am, real & ac, real & ap
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
  real dxm1, dxp1, dxm2, dxp2, lf;      // used for micro resion 
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
      fdm = maxr(fdm,epsl);
      if((clm-0.5)*(clc-0.5)>=0){
        dxm = dxm * fdm;
      } else {
        dxm1 = dxm * fdm;
	aflagm = 0.0;
        tm = tsat;
      }
      if (aflagm==0) {
	if (clm<0.5) {
          lf = lambdal;
          dxm2 = clm * 2.0*(dxm-dxm1);
	} else {
	  lf = lambdav;
          dxm2 = (1.0-clm) * 2.0*(dxm-dxm1);
	}
	/* FDM */
        am = lc*vol*2.0*lf/dxm2/(lc/dxm1+lf/dxm2)/dxm1/(dxm1+dxp);
        ac = lc*vol*2.0/(dxm1*dxp)
           - lc*vol*2.0*lc/dxm1/(lc/dxm1+lf/dxm2)/dxm1/(dxm1+dxp);
        ap = lc*vol*2.0/(dxp*(dxm1+dxp));
        //std::cout<<"f-s-s-FDM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
      } else { 
	/* FVM */
        am = lc * area / dxm * fdm*lm/((1.0-fdm)*lc+fdm*lm);
        ap = 0.5 * (lc + lp) * area / dxp;
        ac = am + ap;
        //std::cout<<"f-s-s-FVM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
      }
    } else if(ofm && onp){
      // s-s-f
      fdp = maxr(fdp,epsl);
      if((clc-0.5)*(clp-0.5)>=0){
        dxp = dxp * fdp;
      } else {
	dxp1 = dxp * fdp;
        aflagp=0.0;
        tp = tsat;
      }
      if (aflagp==0){
	if (clp<0.5) {
          lf = lambdal;
          dxp2 = clp * 2.0*(dxp-dxp1);
	} else {
          lf = lambdav;
          dxp2 = (1.0-clp) * 2.0*(dxp-dxp1);
        }
	real R = dxp2/lf;
	if (Ri>R) {
	  std::cout<<"R:dz/lambdaf= "<<R<<"\n";
          std::cout<<"diff_matrix: need to be develop!!!\n";
          std::cout<<"solid-solid-liquid-interface.\n";
	  exit(0);
	}
        /* FDM */
        am = lc*vol*2.0/(dxm*(dxm+dxp1));
        ac = lc*vol*2.0/(dxm*dxp1)
           - lc*vol*2.0*lc/dxp1/(lf/dxp2+lc/dxp1)/dxp1/(dxm+dxp1);
        ap = lc*vol*2.0*lf/dxp2/(lf/dxp2+lc/dxp1)/dxp1/(dxm+dxp1);
        //std::cout<<"s-s-f-FDM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
      } else {
	/* FVM */
      am = 0.5 * (lc + lm) * area / dxm;
      ap = lc * area / dxp * fdp*lp/((1.0-fdp)*lc+fdp*lp);
      ac = am + ap;
      //std::cout<<"s-s-f-FVM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
      }
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
        dxm=maxr((0.5-clc)/(clm-clc),epsl)*dxm;
        aflagm=0.0;
        tm = tsat;
      }
      if((clc-0.5)*(clp-0.5)>=0){
        dxp=dxp;
      } else {
        dxp=maxr((0.5-clc)/(clp-clc),epsl)*dxp;
        aflagp=0.0;
        tp = tsat;
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
      if ((clc-0.5)*(clm-0.5)>=0) {
        fdm = maxr(fdm,epsl);
        dxm = dxm * fdm;
        if((clc-0.5)*(clp-0.5)>=0){
          dxp=dxp;
        } else {
          dxp=maxr((0.5-clc)/(clp-clc),epsl)*dxp;
          aflagp=0.0;
          tp = tsat;
        }
        if(aflagp==0.0) {
          /* FDM */
          am = lc*vol*2.0/(dxm*(dxm+dxp))*fdm*lm/(fdm*lm+(1.0-fdm)*lc);
          ac = lc*vol*2.0/(dxm*dxp)
             - lc*vol*2.0/(dxm*(dxm+dxp))*(1.0-fdm)*lc/(fdm*lm+(1.0-fdm)*lc);
          ap = lc*vol*2.0/(dxp*(dxm+dxp));
          //std::cout<<"s-f-f-FDM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
        } else {
          /* FVM */
          am = lc * area / dxm * fdm * lm / (fdm*lm+(1.0-fdm)*lc);
          ap = 0.5 * (lc + lp) * area / dxp;
          ac = am + ap;
          //std::cout<<"s-f-f-FVM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
        }
      } else {
	aflagm=0.0;
	if (clc<0.5){
	  dxm = clc*2.0*(dxm-dxm*fdm);
        } else{
	  dxm = (1.0-clc)*2.0*(dxm-dxm*fdm);
        }
	tm = tsat;
        if ((clc-0.5)*(clp-0.5)<0) {
	  dxp = maxr((0.5-clc)/(clp-clc),epsl) * dxp;
	  aflagp = 0.0;
	  tp = tsat;
	}
        am = lc*vol*2.0/dxm/(dxm+dxp);
	ac = lc*vol*2.0/dxm/dxp;
        ap = lc*vol*2.0/dxp/(dxm+dxp);	
      }
    } else if(onm && ofp){
      if ((0.5-clc)/(clp-clc)>=0) {
        // f-f-s
        fdp = maxr(fdp,epsl);
        dxp = dxp * fdp;
        if((clm-0.5)*(clc-0.5)>=0){
          dxm=dxm;
        } else {
          dxm=maxr((0.5-clc)/(clm-clc),epsl)*dxm;
          aflagm=0.0;
          tm = tsat;
        }
        if (aflagm==0.0) {
          /* FDM */
          am = lc*vol*2.0/(dxm*(dxm+dxp));
          ac = lc*vol*2.0/(dxm*dxp)
             - lc*vol*2.0/(dxp*(dxm+dxp))*(1.0-fdp)*lc/((1.0-fdp)*lc+fdp*lp);
          ap = lc*vol*2.0/(dxp*(dxm+dxp))*fdp*lp/((1.0-fdp)*lc+fdp*lp);
          //std::cout<<"f-f-s-FDM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
        } else {
          /* FVM */
          am = 0.5 * (lc + lm) * area / dxm;
          ap = lc * area / dxp * fdp * lp / (fdp*lp+(1.0-fdp)*lc);
          ac = am + ap;
          //std::cout<<"f-f-s-FVM: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
        }
        //std::cout<<"f-f-s: "<<i<<" "<<j<<" "<<k<<" "<<am-ac+ap<<"\n";
      } else {
        aflagp=0.0;
	if (clc<0.5){
	  dxp = clc*2.0*(dxp-dxp*fdp);
        } else{
	  dxp = (1.0-clc)*2.0*(dxp-dxp*fdp);
        }
        tp = tsat;
        if ((clc-0.5)*(clm-0.5)<0) {
          dxm = maxr((0.5-clc)/(clm-clc),epsl) * dxm;
          aflagm = 0.0;
          tm = tsat;
        }
        am = lc*vol*2.0/dxm/(dxm+dxp);
        ac = lc*vol*2.0/dxm/dxp;
        ap = lc*vol*2.0/dxp/(dxm+dxp);
      }
    } else {

      // s-f-s
      std::cout<<"diff_matrix: need to be develop!!!\n";
      std::cout<<"solid-fluid-solid.\n";
      exit(0);

    }
  }

  return;
}
