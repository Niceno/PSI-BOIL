#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Evaluate diffusion for one fluid cell.
*******************************************************************************/
void EnthalpyFD::cell_diffusion_fluid(const Comp m,
                          const int i, const int j, const int k,
                          const int ox, const int oy, const int oz,
                          const real x0, 
                          const coef_gen coef_m, const coef_gen coef_p,
                          const ResistEval re, const Old old, 
                          const real tscn, const real tscm,
                          std::vector<StencilPoint> & stencil,
                          real & Aw, real & Ac, real & Ae, real & F,
                          const Scalar * diff_eddy) {

  stencil.clear();
  std::array<ConnectType,3> ctype = { ConnectType::fluid,
                                      ConnectType::fluid,
                                      ConnectType::fluid };

  real lc = cht.lambda(i,j,k,diff_eddy);
  const real vol = phi.dV(i,j,k);
  real xm,xp,pm,pp;
  Sign markerm, markerp;
  real resinvm(boil::unreal), resinvp(boil::unreal);

  /* is there an interface in west? */
  if(!cht.interface(Sign::neg(),m,i,j,k)) {
    xm=cht.distance_center(Sign::neg(),m,i,j,k);
    pm=phi[i-ox][j-oy][k-oz];
    ctype[0] = ConnectType::fluid;
  } else {
    xm = cht.distance_int(Sign::neg(),m,i,j,k,pm,markerm,re,old);
    ctype[0] = ConnectType::interface;
    /* interfacial resistance plays a role only for liquid */
    if(  cht.use_int_resistance()&&old==Old::no
       &&cht.topo->above_interface(i,j,k)) {
      /* evaluate resistance in the correct cell */
      resinvm = markerm < 0 ?
                cht.evaluate_resinv(m,i   ,j   ,k) :
                cht.evaluate_resinv(m,i-ox,j-oy,k-oz);
    }
  }

  /* is there an interface in east? */
  if(!cht.interface(Sign::pos(),m,i,j,k)) {
    xp=cht.distance_center(Sign::pos(),m,i,j,k);
    pp=phi[i+ox][j+oy][k+oz];
    ctype[2] = ConnectType::fluid;
  } else {
    xp = cht.distance_int(Sign::pos(),m,i,j,k,pp,markerp,re,old);
    ctype[2] = ConnectType::interface;
    /* interfacial resistance plays a role only for liquid */
    if(  cht.use_int_resistance()&&old==Old::no
       &&cht.topo->above_interface(i,j,k)) {
      /* evaluate resistance in the correct cell */
      resinvp = markerp < 0 ?
                cht.evaluate_resinv(m,i   ,j   ,k) :
                cht.evaluate_resinv(m,i+ox,j+oy,k+oz);
    }
  }

  /* matrix coefficients */
  real cxm = (this->*coef_m)(xm,xp,x0) * lc * vol;
  real cxp = (this->*coef_p)(xm,xp,x0) * lc * vol;

  if(old==Old::no) {
    /* build stencil for flux evaluation */
    stencil.push_back(StencilPoint(0,pm,-xm));
    stencil.push_back(StencilPoint(1,phi[i][j][k],0.));
    stencil.push_back(StencilPoint(2,pp, xp));

    /* construct matrix */
    kernel_fluid(ctype,tscn*cxm,tscn*cxp,
                 stencil,resinvm,resinvp,
                 Aw,Ac,Ae,F);
  } else {
    /* simply evaluate */
    Aw = tscm * cxm;
    Ac = tscm * (cxm+cxp);
    Ae = tscm * cxp;
    F += Aw*pm - Ac*phi[i][j][k] + Ae*pp;
  }

  return;
}
