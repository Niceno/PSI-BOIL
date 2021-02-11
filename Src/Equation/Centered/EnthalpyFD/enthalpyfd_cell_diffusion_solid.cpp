#include "enthalpyfd.h"

#include "def.h"

/***************************************************************************//**
*  \brief Evaluate diffusion for one solid cell.
*******************************************************************************/
void EnthalpyFD::cell_diffusion_solid(const Comp m,
                          const int i, const int j, const int k,
                          const int ox, const int oy, const int oz,
                          const real x0, 
                          const coef_gen coef_m, const coef_gen coef_p,
                          const real dSm, const real dSp,
                          const ResistEval re, const Old old, 
                          const real tscn, const real tscm,
                          std::vector<StencilPoint> & stencil,
                          real & Aw, real & Ac, real & Ae, real & F,
                          const Scalar * diff_eddy) {

  Sign dummy; /* dummy sign */
  stencil.clear();
  std::array<ConnectType,3> ctype = { ConnectType::solid,
                                      ConnectType::solid,
                                      ConnectType::solid };

  real lc = cht.lambda(i,j,k,diff_eddy);
  real lcm, lcp;
  const real vol = phi.dV(i,j,k);
  real xm,xp,pm,pp;
  real old_pm,old_pp;
  real dwsrcm, dwsrcp;

  std::array<real,3> resistvals = { boil::unreal,
                                    /* distance to both faces is identical */
                                    cht.distance_face(Sign::neg(),m,i,j,k)/lc,
                                    boil::unreal };

  /* solid also in west? */
  if(dom->ibody().off(i-ox,j-oy,k-oz)) {
    xm=cht.distance_center(Sign::neg(),m,i,j,k);
    old_pm=pm=phi[i-ox][j-oy][k-oz];
    lcm = 0.5*(lc+cht.lambda(i-ox,j-oy,k-oz,diff_eddy));
    ctype[0] = ConnectType::solid;
  } else {
    xm = cht.distance_face(Sign::neg(),m,i,j,k);
    lcm = lc;
    old_pm = cht.node_tmp_sol()[m][i][j][k];
    dwsrcm = cht.dirac_wall_source(i-ox,j-oy,k-oz);

    /* is there an interface? */
    if(!cht.interface(Sign::neg(),m,i,j,k)) {
      real dxfull = cht.distance_center(Sign::neg(),m,i,j,k);
      pm=phi[i-ox][j-oy][k-oz];
      ctype[0] = ConnectType::fluid;

      /* evaluate fluid side resistance */
      resistvals[0] = (dxfull-xm)/cht.lambda(i-ox,j-oy,k-oz,diff_eddy)
                    + cht.wall_resistance(i-ox,j-oy,k-oz);

    } else {
      /* evaluate distance to interface */
      real dxfull = cht.distance_int(Sign::neg(),m,i,j,k,pm,dummy,
                                     ResistEval::no,old);
      ctype[0] = ConnectType::interface;
 
      /* evaluate fluid side resistance, note the inversion of lambda */
      resistvals[0] = (dxfull-xm)/cht.lambda_inv(i-ox,j-oy,k-oz,diff_eddy)
                    + cht.wall_resistance(i-ox,j-oy,k-oz);

      /* interfacial resistance addition */
      if(  cht.use_int_resistance()
         &&cht.topo->below_interface(i-ox,j-oy,k-oz)) {
         resistvals[0] += cht.int_resistance_liq(i-ox,j-oy,k-oz); 
       }
    } 
  }

  /* solid also in east? */
  if(dom->ibody().off(i+ox,j+oy,k+oz)) {
    xp=cht.distance_center(Sign::pos(),m,i,j,k);
    old_pp=pp=phi[i+ox][j+oy][k+oz];
    lcp = 0.5*(lc+cht.lambda(i+ox,j+oy,k+oz,diff_eddy));
    ctype[2] = ConnectType::solid;
  } else {
    xp = cht.distance_face(Sign::pos(),m,i,j,k);
    lcp = lc;
    old_pp = cht.node_tmp_sol()[m][i+ox][j+oy][k+oz];
    dwsrcp = cht.dirac_wall_source(i+ox,j+oy,k+oz);

    /* is there an interface? */
    if(!cht.interface(Sign::pos(),m,i,j,k)) {
      real dxfull = cht.distance_center(Sign::pos(),m,i,j,k);
      pp=phi[i+ox][j+oy][k+oz];
      ctype[2] = ConnectType::fluid;

      /* evaluate fluid side resistance */
      resistvals[2] = (dxfull-xp)/cht.lambda(i+ox,j+oy,k+oz,diff_eddy)
                    + cht.wall_resistance(i+ox,j+oy,k+oz);

    } else {
      /* evaluate distance to interface */
      real dxfull = cht.distance_int(Sign::pos(),m,i,j,k,pp,dummy,
                                     ResistEval::no,old);
      ctype[2] = ConnectType::interface;

      /* evaluate fluid side resistance, note the inversion of lambda */
      resistvals[2] = (dxfull-xp)/cht.lambda_inv(i+ox,j+oy,k+oz,diff_eddy)
                    + cht.wall_resistance(i+ox,j+oy,k+oz);

      /* interfacial resistance addition */
      if(  cht.use_int_resistance()
         &&cht.topo->below_interface(i+ox,j+oy,k+oz)) {
         resistvals[2] += cht.int_resistance_liq(i+ox,j+oy,k+oz);
       }
    }
  }

#ifdef USE_FDM_SOLID
  real cxm = (this->*coef_m)(xm,xp,x0) * lc * vol;
  real cxp = (this->*coef_p)(xm,xp,x0) * lc * vol;
#else
  real cxm = lcm * dSm / xm;
  real cxp = lcp * dSp / xp;
#endif

  if(old==Old::no) {
    kernel_solid(ctype,tscn*cxm,tscn*cxp,pm,pp,
                 resistvals,tscn*dwsrcm,tscn*dwsrcp,
                 Aw,Ac,Ae,F);
  } else {
    Aw = tscm * cxm;
    Ac = tscm * (cxm+cxp);
    Ae = tscm * cxp;
    F += Aw*old_pm - Ac*phi[i][j][k] + Ae*old_pp
       + tscm*(dwsrcm+dwsrcp);
  }

  return;
}
