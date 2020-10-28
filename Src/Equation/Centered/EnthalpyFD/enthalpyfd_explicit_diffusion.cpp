#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Adds diffusion to right hand side.
*
*  This routine adds diffusive terms into right hand side, depending on the
*  time-stepping scheme being used. 
*
*  As an example, for Crank-Nicolson scheme, it is:
*  \f[
*      \{f\}_{old} = \{f\}_{old} + \frac{1}{2} \{D\}^{N-1}
*  \f] 
*  where \f$ \{D\} \f$ is the diffusive term and \f$ N-1 \f$ is the old time 
   step.
*
*  Diffusive term is defined as:
*  \f[
*      \{D\}^{N-1} = [A] \cdot \{ \phi \}^{N-1}
*  \f] 
*  where \f$ [A] \f$ is the system matrix containing only diffusive 
*  contributions and \f$ \{ \phi \}^{N-1} \f$ is the old vector of unknowns. 
*
*  \warning
*  Matrix \f$ [A] \f$ can change in time due to varying pysical properties.
*  To be fully consistent in such a case, this subroutine has to be invoked
*  before forming the new system of governing equations. This sequence is 
*  stipulated in main program. Obviously, these concerns are important for 
*  non-implicit time stepping schemes only.
*******************************************************************************/
void EnthalpyFD::explicit_diffusion(const Scalar * diff_eddy) {

  phi.exchange();

  /* get time stepping coefficient */
  real tscn = diff_ts.N();
  assert( tscn > 0.0 );
  real tsc = diff_ts.Nm1(); /* 0.5 for c.n. 0.0 for fully implicit */

  if( accelerated_no_solid ) {
    /*------------------------+ 
    |  x direction (w and e)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      const real pc = phi[i][j][k];
      real xm,xp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if(!cht.interface(Sign::neg(),Comp::i(),i,j,k)){
        xm=phi.dxw(i);
        pm=phi[i-1][j][k];
      } else {
        xm = cht.distance_int_x(Sign::neg(),i,j,k,pm);
        aflagm=1.0;
      }
      if(!cht.interface(Sign::pos(),Comp::i(),i,j,k)){
        xp=phi.dxe(i);
        pp=phi[i+1][j][k];
      } else {
        xp = cht.distance_int_x(Sign::pos(),i,j,k,pp);
        aflagp=1.0;
      }
      real cxm = coef_x_m(xm,xp,phi.xc(i));
      real cxp = coef_x_p(xm,xp,phi.xc(i));

      const real Aw = tsc * lc * vol * cxm;
      const real Ac = tsc * lc * vol * (cxm+cxp);
      const real Ae = tsc * lc * vol * cxp;
      fold[i][j][k] += Aw*pm - Ac*pc + Ae*pp;
      /* add implicit part of the diffusion term */
      const real Awi = tscn * lc * vol * cxm * aflagm;
      const real Aei = tscn * lc * vol * cxp * aflagp;
      ftif[i][j][k] = Awi*pm + Aei*pp;
    }

    /*------------------------+ 
    |  y direction (s and n)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      const real pc = phi[i][j][k];
      real ym,yp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if(!cht.interface(Sign::neg(),Comp::j(),i,j,k)){
        ym=phi.dys(j);
        pm=phi[i][j-1][k];
      } else {
        ym = cht.distance_int_y(Sign::neg(),i,j,k,pm);
        aflagm=1.0;
      }
      if(!cht.interface(Sign::pos(),Comp::j(),i,j,k)){
        yp=phi.dyn(j);
        pp=phi[i][j+1][k];
      } else {
        yp = cht.distance_int_y(Sign::pos(),i,j,k,pp);
        aflagp=1.0;
      }
      real cym = coef_y_m(ym,yp,phi.yc(j));
      real cyp = coef_y_p(ym,yp,phi.yc(j));

      const real As = tsc * lc * vol * cym;
      const real Ac = tsc * lc * vol * (cym+cyp);
      const real An = tsc * lc * vol * cyp;
      fold[i][j][k] += As*pm - Ac*pc + An*pp;
      const real Asi = tscn * lc * vol * cym * aflagm;
      const real Ani = tscn * lc * vol * cyp * aflagp;
      ftif[i][j][k] += Asi*pm + Ani*pp;
    }

    /*------------------------+ 
    |  z direction (b and t)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      const real pc = phi[i][j][k];
      real zm,zp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if(!cht.interface(Sign::neg(),Comp::k(),i,j,k)){
        zm=phi.dzb(k);
        pm=phi[i][j][k-1];
      } else {
        zm = cht.distance_int_z(Sign::neg(),i,j,k,pm);
        aflagm=1.0;
      }
      if(!cht.interface(Sign::pos(),Comp::k(),i,j,k)){
        zp=phi.dzt(k);
        pp=phi[i][j][k+1];
      } else {
        zp = cht.distance_int_z(Sign::pos(),i,j,k,pp);
        aflagp=1.0;
      }
      real czm = coef_z_m(zm,zp,phi.zc(k));
      real czp = coef_z_p(zm,zp,phi.zc(k));
      
      const real Ab = tsc * lc * vol * czm;
      const real Ac = tsc * lc * vol * (czm+czp);
      const real At = tsc * lc * vol * czp;
      fold[i][j][k] += Ab*pm - Ac*pc + At*pp;
      const real Abi = tscn * lc * vol * czm * aflagm;
      const real Ati = tscn * lc * vol * czp * aflagp;
      ftif[i][j][k] += Abi*pm + Ati*pp;
    }
    // need to add here immersed boundary without solid !!!

  /*---------------------------------------------+
  |  can feature conduction through solid parts  |
  +---------------------------------------------*/
  } else {

    /* if this is used without solid, that is, if accelerated_no_solid == true
       and at the same time solid() == false, it is necessary to prevent
       segfaults when the 'solid' conductivities are referenced. However, since
       they will not be used for anything inside diff_matrix, we can set any
       value to them. For acceleration (to avoid checking a flag for each cell
       and each direction), the matter pointer is set to fluid in the constructor
       instead. Note that if, at one point, existence of ibodies without solid
       conduction is allowed, this of course needs rewriting */

    ftif = 0.;

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
        bool onm, onc, onp, ofm, ofc, ofp; // on & off
        real lsm, lsc, lsp; // lambdas
        real lvm, lvc, lvp; // lambdav
        real llm, llc, llp; // lambdal
        int clm, clc, clp; // topology flag
        real dxm, dxp, fdm, fdp, fdms, fdps;
        real pm, pc, pp;
        real am, ac, ap;
        real tm, tc, tp;
        real aflagm, aflagp;
        real aream, areap;
        real sourceterm(0.0);
        real pos0;
        coef_gen coef_m, coef_p;

        onm=dom->ibody().on (i-ii,j-jj,k-kk);
        onc=dom->ibody().on (i   ,j   ,k   );
        onp=dom->ibody().on (i+ii,j+jj,k+kk);
        ofm=dom->ibody().off(i-ii,j-jj,k-kk);
        ofc=dom->ibody().off(i   ,j   ,k   );
        ofp=dom->ibody().off(i+ii,j+jj,k+kk);
        lsm=safe_solid->lambda (i-ii,j-jj,k-kk);
        lsc=safe_solid->lambda (i   ,j   ,k   );
        lsp=safe_solid->lambda (i+ii,j+jj,k+kk);
        lvm=cht.lambdav(i-ii,j-jj,k-kk,diff_eddy);
        lvc=cht.lambdav(i   ,j   ,k   ,diff_eddy);
        lvp=cht.lambdav(i+ii,j+jj,k+kk,diff_eddy);
        llm=cht.lambdal(i-ii,j-jj,k-kk,diff_eddy);
        llc=cht.lambdal(i   ,j   ,k   ,diff_eddy);
        llp=cht.lambdal(i+ii,j+jj,k+kk,diff_eddy);
        clm=iflag[i-ii][j-jj][k-kk];
        clc=iflag[i   ][j   ][k   ];
        clp=iflag[i+ii][j+jj][k+kk];
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
        diff_matrix(am, ac, ap
                  , tm, tc, tp
                  , aflagm, aflagp
                  , sourceterm
                  , pos0, coef_m, coef_p
                  , vol, aream, areap
                  , onm, onc, onp, ofm, ofc, ofp
                  , lsm, lsc, lsp
                  , lvm, lvc, lvp
                  , llm, llc, llp
                  , clm, clc, clp
                  , dxm, dxp, fdm, fdp, fdms, fdps
                  , pm, pc, pp
                  , i, j, k, m);
  
        fold[i][j][k] += sourceterm;
        if(tsc!=0){
          fold[i][j][k] += tsc * (am*tm*aflagm - ac*tc + ap*tp*aflagp);
        }
        ftif[i][j][k] += tscn * (am*(1.0-aflagm)*tm+ap*(1.0-aflagp)*tp);
      }
    }
  } /* solid conduction */

  return;

}
