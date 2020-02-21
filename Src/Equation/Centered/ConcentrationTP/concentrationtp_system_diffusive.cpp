#include "concentrationtp.h"

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void ConcentrationTP::create_system_diffusive(const Scalar * diff_eddy) {

  /* initialize: get time stepping coefficient */
  real tscn = diff_ts.N();
  real gm, gp, lm, lp;
  assert( tscn > 0.0 );

  /* coefficients in i direction (w and e) */
  for_ijk(i,j,k) {
    Comp mcomp = Comp::u();

    real clrm = surface_color(Sign::neg(),mcomp,i,j,k);
    real clrp = surface_color(Sign::pos(),mcomp,i,j,k);

    real dxw = phi.dxw(i);
    real dxe = phi.dxe(i);
  
    /* diffusion in gas */
    real dcm = dcoef->value(mcomp,i  ,j,k);
    real dcp = dcoef->value(mcomp,i+1,j,k);

    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dcm += ((*diff_eddy)[i-1][j][k] +
              (*diff_eddy)[i  ][j][k] )/(2.0*turbS);
      dcp += ((*diff_eddy)[i  ][j][k] +
              (*diff_eddy)[i+1][j][k] )/(2.0*turbS);
    } 
    real sw = dSx(i,j,k)*(1.0-clrm);
    real se = dSx(i,j,k)*(1.0-clrp);

    gm =  tscn * dcm * sw / dxw;
    gp =  tscn * dcp * se / dxe;

    A.w[i][j][k] = gm;
    A.e[i][j][k] = gp;
    A.c[i][j][k] += gm + gp;
  }

  /* coefficients in j direction (s and n) */
  for_ijk(i,j,k) {
    Comp mcomp = Comp::v();

    real clrm = surface_color(Sign::neg(),mcomp,i,j,k);
    real clrp = surface_color(Sign::pos(),mcomp,i,j,k);

    real dys = phi.dys(j);
    real dyn = phi.dyn(j);
    
    /* diffusion in gas */
    real dcm = dcoef->value(mcomp,i,j  ,k);
    real dcp = dcoef->value(mcomp,i,j+1,k);

    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dcm += ((*diff_eddy)[i][j-1][k] +
              (*diff_eddy)[i][j  ][k] )/(2.0*turbS);
      dcp += ((*diff_eddy)[i][j  ][k] +
              (*diff_eddy)[i][j+1][k] )/(2.0*turbS);
    } 
    real ss = dSy(i,j,k)*(1.0-clrm);
    real sn = dSy(i,j,k)*(1.0-clrp);

    gm = tscn * dcm * ss / dys;
    gp = tscn * dcp * sn / dyn;

    A.s[i][j][k] = gm;
    A.n[i][j][k] = gp;
    A.c[i][j][k] += gm + gp;
  }

  /* coefficients in k direction (b and t) */
  for_ijk(i,j,k) {
    Comp mcomp = Comp::w();

    real clrm = surface_color(Sign::neg(),mcomp,i,j,k);
    real clrp = surface_color(Sign::pos(),mcomp,i,j,k);

    real dzb = phi.dzb(k);
    real dzt = phi.dzt(k);

    /* diffusion in gas */
    real dcm = dcoef->value(mcomp,i,j,k  );
    real dcp = dcoef->value(mcomp,i,j,k+1);

    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dcm += ((*diff_eddy)[i][j][k-1] +
              (*diff_eddy)[i][j][k  ] )/(2.0*turbS);
      dcp += ((*diff_eddy)[i][j][k  ] +
              (*diff_eddy)[i][j][k+1] )/(2.0*turbS);
    } 
    real sb = dSz(i,j,k)*(1.0-clrm);
    real st = dSz(i,j,k)*(1.0-clrp);

    gm = tscn * dcm * sb / dzb;
    gp = tscn * dcp * st / dzt;

    A.b[i][j][k] = gm;
    A.t[i][j][k] = gp;
    A.c[i][j][k] += gm + gp;
  }

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) {
    for(int cc=0; cc<dom->ibody().nccells(); cc++) {
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

      /* w */
      if( dom->ibody().on(i-1,j,k) ) {
        const real fSw  = dom->ibody().fSw(cc);
        const real fdxw = dom->ibody().fdxw(cc);
        if( dom->ibody().on(i,j,k) ) A.w[i]  [j][k] *= (fSw / fdxw);
        else if(fdxw != 1.0)         A.e[i-1][j][k] /= (1.0 - fdxw);
      }

      /* e */
      if( dom->ibody().on(i+1,j,k) ) {
        const real fSe  = dom->ibody().fSe(cc);
        const real fdxe = dom->ibody().fdxe(cc);
        if( dom->ibody().on(i,j,k) ) A.e[i]  [j][k] *= (fSe / fdxe);
        else if(fdxe != 1.0)         A.w[i+1][j][k] /= (1.0 - fdxe);
      }

      /* s */
      if( dom->ibody().on(i,j-1,k) ) {
        const real fSs  = dom->ibody().fSs(cc);
        const real fdys = dom->ibody().fdys(cc);
        if( dom->ibody().on(i,j,k) ) A.s[i][j]  [k] *= (fSs / fdys);
        else if(fdys != 1.0)         A.n[i][j-1][k] /= (1.0 - fdys);
      }

      /* n */
      if( dom->ibody().on(i,j+1,k) ) {
        const real fSn  = dom->ibody().fSn(cc);
        const real fdyn = dom->ibody().fdyn(cc);
        if( dom->ibody().on(i,j,k) ) A.n[i][j]  [k] *= (fSn / fdyn);
        else if(fdyn != 1.0)         A.s[i][j+1][k] /= (1.0 - fdyn);
      }

      /* b */
      if( dom->ibody().on(i,j,k-1) ) {
        const real fSb  = dom->ibody().fSb(cc);
        const real fdzb = dom->ibody().fdzb(cc);
        if( dom->ibody().on(i,j,k) ) A.b[i][j][k]   *= (fSb / fdzb);
        else if(fdzb != 1.0)         A.t[i][j][k-1] /= (1.0 - fdzb);
      }
      /* t */
      if( dom->ibody().on(i,j,k+1) ) {
        const real fSt  = dom->ibody().fSt(cc);
        const real fdzt = dom->ibody().fdzt(cc);
        if( dom->ibody().on(i,j,k) ) A.t[i][j][k]   *= (fSt / fdzt);
        else if(fdzt != 1.0)         A.b[i][j][k+1] /= (1.0 - fdzt);
      }

    }

  } /* is there an immersed body */
}
