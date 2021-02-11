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

  /* clear the stencil first */
  stencil.clear();

  std::array<ConnectType,3> ctype = { ConnectType::fluid,
                                      ConnectType::fluid,
                                      ConnectType::fluid };

  real lc = cht.lambda(i,j,k,diff_eddy);
  const real vol = phi.dV(i,j,k);
  real xm(boil::unreal),xp(boil::unreal);
  real new_pm(boil::unreal),new_pp(boil::unreal);
  real old_pm(boil::unreal),old_pp(boil::unreal);
  Sign markerm, markerp;
  real resinvm(boil::unreal), resinvp(boil::unreal);
  real dwsrcm(boil::unreal), dwsrcp(boil::unreal);

  std::array<real,3> resistvals = { boil::unreal,
                                    /* distance to both faces is identical */
                                    cht.distance_face(Sign::neg(),m,i,j,k)/lc,
                                    boil::unreal };

  /* is there no interface in west? */
  if(!cht.interface(Sign::neg(),m,i,j,k,old)) {
    xm = cht.distance_center(Sign::neg(),m,i,j,k);
    old_pm = new_pm = phi[i-ox][j-oy][k-oz];
    ctype[0] = ConnectType::fluid;
  /* is there a wall in the west? */
  } else if(dom->ibody().off(i-ox,j-oy,k-oz)) {
    xm = cht.distance_face(Sign::neg(),m,i,j,k);
    new_pm = phi[i-ox][j-oy][k-oz];
    old_pm = cht.node_tmp_flu()[m][i][j][k];
    dwsrcm = cht.dirac_wall_source(i,j,k);
    ctype[0] = ConnectType::solid;

    /* evaluate solid side resistance */
    real dxfull = cht.distance_center(Sign::neg(),m,i,j,k);
    resistvals[0] = (dxfull-xm)/cht.lambda(i-ox,j-oy,k-oz,diff_eddy)
                  + cht.wall_resistance(i,j,k);
  } else {
    xm = cht.distance_int(Sign::neg(),m,i,j,k,new_pm,markerm,re,old);
    old_pm = new_pm;
    ctype[0] = ConnectType::interface;
    /* interfacial resistance plays a role only for liquid */
    if(  cht.use_int_resistance()&&old==Old::no
       &&cht.above_interface(i,j,k,old)) {
      /* evaluate resistance in the correct cell */
      resinvm = markerm < 0 ?
                cht.evaluate_resinv(m,i   ,j   ,k) :
                cht.evaluate_resinv(m,i-ox,j-oy,k-oz);
    }
  }

  /* is there no interface in east? */
  if(!cht.interface(Sign::pos(),m,i,j,k,old)) {
    xp = cht.distance_center(Sign::pos(),m,i,j,k);
    old_pp = new_pp = phi[i+ox][j+oy][k+oz];
    ctype[2] = ConnectType::fluid;
  /* is there a wall in the east */
  } else if(dom->ibody().off(i+ox,j+oy,k+oz)) {
    xp = cht.distance_face(Sign::pos(),m,i,j,k);
    new_pp = phi[i+ox][j+oy][k+oz];
    old_pp = cht.node_tmp_flu()[m][i+ox][j+oy][k+oz];
    dwsrcp = cht.dirac_wall_source(i,j,k);
    ctype[2] = ConnectType::solid;

    /* evaluate solid side resistance */
    real dxfull = cht.distance_center(Sign::pos(),m,i,j,k);
    resistvals[2] = (dxfull-xp)/cht.lambda(i+ox,j+oy,k+oz,diff_eddy)
                    + cht.wall_resistance(i,j,k);
  } else {
    xp = cht.distance_int(Sign::pos(),m,i,j,k,new_pp,markerp,re,old);
    old_pp = new_pp;
    ctype[2] = ConnectType::interface;
    /* interfacial resistance plays a role only for liquid */
    if(  cht.use_int_resistance()&&old==Old::no
       &&cht.above_interface(i,j,k,old)) {
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
    stencil.push_back(StencilPoint(0,new_pm,-xm));
    stencil.push_back(StencilPoint(1,phi[i][j][k],0.));
    stencil.push_back(StencilPoint(2,new_pp, xp));

    /* construct matrix -- accelerated kernel without solid */
    if(ctype[0] != ConnectType::solid && ctype[2] != ConnectType::solid) {
      kernel_fluid1(ctype,tscn*cxm,tscn*cxp,
                    stencil,resinvm,resinvp,
                    Aw,Ac,Ae,F);
    /* kernel with solid */
    } else {
      kernel_fluid2(ctype,tscn*cxm,tscn*cxp,
                    stencil,resinvm,resinvp,
                    resistvals,dwsrcm,dwsrcp,
                    Aw,Ac,Ae,F);
    }
  } else {
    /* simply evaluate */
    Aw = tscm * cxm;
    Ac = tscm * (cxm+cxp);
    Ae = tscm * cxp;
    F += Aw*old_pm - Ac*phi[i][j][k] + Ae*old_pp;
  }

  return;
}
