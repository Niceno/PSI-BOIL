#include "enthalpyfd.h"

using namespace std;

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
void EnthalpyFD::diffusion_fd(const Scalar * diff_eddy) {

  phi.exchange();

  /* get time stepping coefficient */
  real tscn = diff_ts.N();
  assert( tscn > 0.0 );
  real tsc = diff_ts.Nm1(); /* 0.5 for c.n. 0.0 for fully implicit */

  if( !solid() ) {
    /*------------------------+ 
    |  x direction (w and e)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc, cp_mass;
      if((*clr)[i][j][k]>=0.5){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrw = std::min(std::max((*clr)[i-1][j][k],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i  ][j][k],0.0),1.0);
      const real clre = std::min(std::max((*clr)[i+1][j][k],0.0),1.0);
      const real pc = phi[i][j][k];
      real xm,xp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if((clrw-0.5)*(clrc-0.5)>=0){
        xm=phi.dxw(i);
        pm=phi[i-1][j][k];
      } else {
        xm=std::max((0.5-clrc)/(clrw-clrc),epsl)*phi.dxw(i);
        pm=tsat;
        aflagm=1.0;
      }
      if((clrc-0.5)*(clre-0.5)>=0){
        xp=phi.dxe(i);
        pp=phi[i+1][j][k];
      } else {
        xp=std::max((0.5-clrc)/(clre-clrc),epsl)*phi.dxe(i);
        pp=tsat;
        aflagp=1.0;
      }
      const real Aw = tsc * lc * vol * 2.0 / (xm*(xm+xp));
      const real Ac = tsc * lc * vol * 2.0 / (xm*xp);
      const real Ae = tsc * lc * vol * 2.0 / (xp*(xm+xp));
      fold[i][j][k] += Aw*pm - Ac*pc + Ae*pp;
      /* add implicit part of the diffusion term */
      const real Awi = tscn * lc * vol * 2.0 / (xm*(xm+xp)) * aflagm;
      const real Aei = tscn * lc * vol * 2.0 / (xp*(xm+xp)) * aflagp;
      fold[i][j][k] += Awi*pm + Aei*pp;
    }

    /*------------------------+ 
    |  y direction (s and n)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc,cp_mass;
      if((*clr)[i][j][k]>=0.5){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrs = std::min(std::max((*clr)[i][j-1][k],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j  ][k],0.0),1.0);
      const real clrn = std::min(std::max((*clr)[i][j+1][k],0.0),1.0);
      const real pc = phi[i][j][k];
      real ym,yp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if((clrs-0.5)*(clrc-0.5)>=0){
        ym=phi.dys(j);
        pm=phi[i][j-1][k];
      } else {
        ym=std::max((0.5-clrc)/(clrs-clrc),epsl)*phi.dys(j);
        pm=tsat;
        aflagm=1.0;
      }
      if((clrc-0.5)*(clrn-0.5)>=0){
        yp=phi.dyn(j);
        pp=phi[i][j+1][k];
      } else {
        yp=std::max((0.5-clrc)/(clrn-clrc),epsl)*phi.dyn(j);
        pp=tsat;
        aflagp=1.0;
      }
      const real As = tsc * lc * vol * 2.0 / (ym*(ym+yp));
      const real Ac = tsc * lc * vol * 2.0 / (ym*yp);
      const real An = tsc * lc * vol * 2.0 / (yp*(ym+yp));
      fold[i][j][k] += As*pm - Ac*pc + An*pp;
      const real Asi = tscn * lc * vol * 2.0 / (ym*(ym+yp)) * aflagm;
      const real Ani = tscn * lc * vol * 2.0 / (yp*(ym+yp)) * aflagp;
      fold[i][j][k] += Asi*pm + Ani*pp;
    }

    /*------------------------+ 
    |  z direction (b and t)  |
    +------------------------*/
    for_ijk(i,j,k){
      real lc,cp_mass;
      if((*clr)[i][j][k]>=0.5){
        lc = lambdal;
        cp_mass = cpl/rhol;
      } else {
        lc = lambdav;
        cp_mass = cpv/rhov;
      }
      if(diff_eddy){
        lc += (*diff_eddy)[i][j][k]*cp_mass/turbP;
      }
      const real vol = phi.dV(i,j,k);
      const real clrb = std::min(std::max((*clr)[i][j][k-1],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j][k  ],0.0),1.0);
      const real clrt = std::min(std::max((*clr)[i][j][k+1],0.0),1.0);
      const real pc = phi[i][j][k];
      real zm,zp,pm,pp,aflagm,aflagp;
      aflagm=aflagp=0.0;
      if((clrb-0.5)*(clrc-0.5)>=0){
        zm=phi.dzb(k);
        pm=phi[i][j][k-1];
      } else {
        zm=std::max((0.5-clrc)/(clrb-clrc),epsl)*phi.dzb(k);
        pm=tsat;
        aflagm=1.0;
      }
      if((clrc-0.5)*(clrt-0.5)>=0){
        zp=phi.dzt(k);
        pp=phi[i][j][k+1];
      } else {
        zp=std::max((0.5-clrc)/(clrt-clrc),epsl)*phi.dzt(k);
        pp=tsat;
        aflagp=1.0;
      }
      const real Ab = tsc * lc * vol * 2.0 / (zm*(zm+zp));
      const real Ac = tsc * lc * vol * 2.0 / (zm*zp);
      const real At = tsc * lc * vol * 2.0 / (zp*(zm+zp));
      fold[i][j][k] += Ab*pm - Ac*pc + At*pp;
      const real Abi = tscn * lc * vol * 2.0 / (zm*(zm+zp)) * aflagm;
      const real Ati = tscn * lc * vol * 2.0 / (zp*(zm+zp)) * aflagp;
      fold[i][j][k] += Abi*pm + Ati*pp;
    }
    // need to add here immersed boundary without solid !!!

#if 1
  /*------------------------------------------+
  |  features conduction through solid parts  |
  +------------------------------------------*/
  } else {

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
        real lc, xm, xp;

        bool onm, onc, onp, ofm, ofc, ofp; // on & off
        real lsm, lsc, lsp, lfm, lfc, lfp; // lambda
        real clm, clc, clp; // color function
        real dxm, dxp, fdm, fdp, fdms, fdps;
        real edm, edc, edp; // eddy viscosity
        real pm, pc, pp;
        real am, ac, ap;
        real tm, tc, tp;
        real aflagm, aflagp;
        real area;

        onm=dom->ibody().on (i-ii,j-jj,k-kk);
        onc=dom->ibody().on (i   ,j   ,k   );
        onp=dom->ibody().on (i+ii,j+jj,k+kk);
        ofm=dom->ibody().off(i-ii,j-jj,k-kk);
        ofc=dom->ibody().off(i   ,j   ,k   );
        ofp=dom->ibody().off(i+ii,j+jj,k+kk);
        lsm=solid()->lambda (i-ii,j-jj,k-kk);
        lsc=solid()->lambda (i   ,j   ,k   );
        lsp=solid()->lambda (i+ii,j+jj,k+kk);
        lfm=fluid()->lambda (i-ii,j-jj,k-kk);
        lfc=fluid()->lambda (i   ,j   ,k   );
        lfp=fluid()->lambda (i+ii,j+jj,k+kk);
        clm=(*clr)[i-ii][j-jj][k-kk];
        clc=(*clr)[i   ][j   ][k   ];
        clp=(*clr)[i+ii][j+jj][k+kk];
	if(m==Comp::i()){
          dxm=phi.dxw(i);
          dxp=phi.dxe(i);
          area = dSx(i,j,k);
          fdm=dom->ibody().fdxw(i,j,k);
          fdp=dom->ibody().fdxe(i,j,k);
          fdms=dom->ibody().fdxe(i-1,j,k);
          fdps=dom->ibody().fdxw(i+1,j,k);
        } else if(m==Comp::j()){
          dxm=phi.dys(j);
          dxp=phi.dyn(j);
          area = dSy(i,j,k);
          fdm=dom->ibody().fdys(i,j,k);
          fdp=dom->ibody().fdyn(i,j,k);
          fdms=dom->ibody().fdyn(i,j-1,k);
          fdps=dom->ibody().fdys(i,j+1,k);
        } else {
          dxm=phi.dzb(k);
          dxp=phi.dzt(k);
          area = dSz(i,j,k);
          fdm=dom->ibody().fdzb(i,j,k);
          fdp=dom->ibody().fdzt(i,j,k);
          fdms=dom->ibody().fdzt(i,j,k-1);
          fdps=dom->ibody().fdzb(i,j,k+1);
        }
        pm=phi[i-ii][j-jj][k-kk];
        pc=phi[i   ][j   ][k   ];
        pp=phi[i+ii][j+jj][k+kk];
        edm=edc=edp=0.0;
        if(diff_eddy){
          if (onm) edm = (*diff_eddy)[i-ii][j-jj][k-kk];
          if (onc) edc = (*diff_eddy)[i   ][j   ][k   ];
          if (onp) edp = (*diff_eddy)[i+ii][j+jj][k+kk];
        }
        diff_matrix(am, ac, ap
                  , tm, tc, tp
                  , aflagm, aflagp
                  , vol, area
                  , onm, onc, onp, ofm, ofc, ofp
                  , lsm, lsc, lsp, lfm, lfc, lfp
                  , clm, clc, clp
                  , dxm, dxp, fdm, fdp, fdms, fdps
                  , pm, pc, pp
                  , edm, edc, edp
                  , i, j, k, m);
  
        if(tsc!=0){
          fold[i][j][k] += tsc * (am*tm*aflagm - ac*tc + ap*tp*aflagp);
        }
        fold[i][j][k] += tscn * (am*(1.0-aflagm)*tm+ap*(1.0-aflagp)*tp);
      }
    }
  }
#endif
}
/*-----------------------------------------------------------------------------+
 '$Id: enthalpyfd_diffusion_fd.cpp,v 1.6 2015/06/29 18:21:19 sato Exp $'/
+-----------------------------------------------------------------------------*/
