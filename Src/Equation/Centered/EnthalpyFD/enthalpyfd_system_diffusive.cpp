#include "enthalpyfd.h"
using namespace std;

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void EnthalpyFD::create_system_diffusive(const Scalar * diff_eddy) {

  /* initialize: get time stepping coefficient */
  real tscn = diff_ts.N();
  assert( tscn > 0.0 );

  /*------------------------------------+
  |  no conduction through solid parts  |
  +------------------------------------*/
  if( !solid() ) {

    /* coefficients in i direction (w and e) */
    for_ijk(i,j,k) {
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
      real xm,xp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if((clrw-0.5)*(clrc-0.5)>=0){
        xm=phi.dxw(i);
      } else {
        xm=std::max((0.5-clrc)/(clrw-clrc),epsl)*phi.dxw(i);
        aflagm=0.0;
      }
      if((clrc-0.5)*(clre-0.5)>=0){
        xp=phi.dxe(i);
      } else {
        xp=std::max((0.5-clrc)/(clre-clrc),epsl)*phi.dxe(i);
        aflagp=0.0;
      }
      A.w[i][j][k] =  tscn * lc * vol * 2.0 / (xm*(xm+xp)) * aflagm;
      A.c[i][j][k] += tscn * lc * vol * 2.0 / (xm*xp);
      A.e[i][j][k] =  tscn * lc * vol * 2.0 / (xp*(xm+xp)) * aflagp;
    }

    /* coefficients in j direction (s and n) */
    for_ijk(i,j,k) {
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
      const real clrs = std::min(std::max((*clr)[i][j-1][k],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j  ][k],0.0),1.0);
      const real clrn = std::min(std::max((*clr)[i][j+1][k],0.0),1.0);
      real ym,yp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if((clrs-0.5)*(clrc-0.5)>=0){
        ym=phi.dys(j);
      } else {
        ym=std::max((0.5-clrc)/(clrs-clrc),epsl)*phi.dys(j);
        aflagm=0.0;
      }
      if((clrc-0.5)*(clrn-0.5)>=0){
        yp=phi.dyn(j);
      } else {
        yp=std::max((0.5-clrc)/(clrn-clrc),epsl)*phi.dyn(j);
        aflagp=0.0;
      }
      A.s[i][j][k] =  tscn * lc * vol * 2.0 / (ym*(ym+yp)) * aflagm;
      A.c[i][j][k] += tscn * lc * vol * 2.0 / (ym*yp);
      A.n[i][j][k] =  tscn * lc * vol * 2.0 / (yp*(ym+yp)) * aflagp;
    }

    /* coefficients in k direction (b and t) */
    for_ijk(i,j,k) {
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
      const real clrb = std::min(std::max((*clr)[i][j][k-1],0.0),1.0);
      const real clrc = std::min(std::max((*clr)[i][j][k  ],0.0),1.0);
      const real clrt = std::min(std::max((*clr)[i][j][k+1],0.0),1.0);
      real zm,zp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if((clrb-0.5)*(clrc-0.5)>=0){
        zm=phi.dzb(k);
      } else {
        zm=std::max((0.5-clrc)/(clrb-clrc),epsl)*phi.dzb(k);
        aflagm=0.0;
      }
      if((clrc-0.5)*(clrt-0.5)>=0){
        zp=phi.dzt(k);
      } else {
        zp=std::max((0.5-clrc)/(clrt-clrc),epsl)*phi.dzt(k);
        aflagp=0.0;
      }
      A.b[i][j][k] =  tscn * lc * vol * 2.0 / (zm*(zm+zp)) * aflagm;
      A.c[i][j][k] += tscn * lc * vol * 2.0 / (zm*zp);
      A.t[i][j][k] =  tscn * lc * vol * 2.0 / (zp*(zm+zp)) * aflagp;
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

  /*------------------------------------------+
  |  features conduction through solid parts  |
  +------------------------------------------*/
  } else {

    for_ijk(i,j,k) {

      const real vol = phi.dV(i,j,k);
      bool onm, onc, onp, ofm, ofc, ofp; // on & off
      real lsm, lsc, lsp, lfm, lfc, lfp; // lambda
      real clm, clc, clp; // color function
      real dxm, dxp, fdm, fdp, fdms, fdps;
      real pm, pc, pp;    // temperature-input
      real edm, edc, edp; // eddy viscosity
      real am, ac, ap;
      real tm, tc, tp;    // temperature-output
      real aflagm, aflagp;

      /* i-direction */
      onm=dom->ibody().on(i-1,j,k);
      onc=dom->ibody().on(i  ,j,k);
      onp=dom->ibody().on(i+1,j,k);
      ofm=dom->ibody().off(i-1,j,k);
      ofc=dom->ibody().off(i  ,j,k);
      ofp=dom->ibody().off(i+1,j,k);
      lsm=solid()->lambda(i-1,j,k);
      lsc=solid()->lambda(i  ,j,k);
      lsp=solid()->lambda(i+1,j,k);
      lfm=fluid()->lambda(i-1,j,k);
      lfc=fluid()->lambda(i  ,j,k);
      lfp=fluid()->lambda(i+1,j,k);
      clm=(*clr)[i-1][j][k];
      clc=(*clr)[i  ][j][k];
      clp=(*clr)[i+1][j][k];
      dxm=phi.dxw(i);
      dxp=phi.dxe(i);
      fdm=dom->ibody().fdxw(i,j,k);
      fdp=dom->ibody().fdxe(i,j,k);
      fdms=dom->ibody().fdxe(i-1,j,k);
      fdps=dom->ibody().fdxw(i+1,j,k);
      pm=phi[i-1][j][k];
      pc=phi[i  ][j][k];
      pp=phi[i+1][j][k];
      edm=edc=edp=0.0;
      if(diff_eddy){
        if (onm) edm = (*diff_eddy)[i-1][j][k];
        if (onc) edc = (*diff_eddy)[i  ][j][k];
        if (onp) edp = (*diff_eddy)[i+1][j][k];
      }

      diff_matrix(am, ac, ap
                , tm, tc, tp
                , aflagm, aflagp
                , vol, dSx(i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp, lfm, lfc, lfp
                , clm, clc, clp
                , dxm, dxp, fdm, fdp, fdms, fdps
                , pm, pc, pp
                , edm, edc, edp
                , i, j, k, Comp::u());
      A.w[i][j][k] = tscn * am * aflagm;
      A.c[i][j][k]+= tscn * ac;
      A.e[i][j][k] = tscn * ap * aflagp;

      /* j-direction */
      onm=dom->ibody().on(i,j-1,k);
      onc=dom->ibody().on(i,j  ,k);
      onp=dom->ibody().on(i,j+1,k);
      ofm=dom->ibody().off(i,j-1,k);
      ofc=dom->ibody().off(i,j  ,k);
      ofp=dom->ibody().off(i,j+1,k);
      lsm=solid()->lambda(i,j-1,k);
      lsc=solid()->lambda(i,j  ,k);
      lsp=solid()->lambda(i,j+1,k);
      lfm=fluid()->lambda(i,j-1,k);
      lfc=fluid()->lambda(i,j  ,k);
      lfp=fluid()->lambda(i,j+1,k);
      clm=(*clr)[i][j-1][k];
      clc=(*clr)[i][j  ][k];
      clp=(*clr)[i][j+1][k];
      dxm=phi.dys(j);
      dxp=phi.dyn(j);
      fdm=dom->ibody().fdys(i,j,k);
      fdp=dom->ibody().fdyn(i,j,k);
      fdms=dom->ibody().fdyn(i,j-1,k);
      fdps=dom->ibody().fdys(i,j+1,k);
      pm=phi[i][j-1][k];
      pc=phi[i][j  ][k];
      pp=phi[i][j+1][k];
      edm=edc=edp=0.0;
      if(diff_eddy){
        if (onm) edm = (*diff_eddy)[i][j-1][k];
        if (onc) edc = (*diff_eddy)[i][j  ][k];
        if (onp) edp = (*diff_eddy)[i][j+1][k];
      }
  
      diff_matrix(am, ac, ap
                , tm, tc, tp
                , aflagm, aflagp
                , vol, dSy(i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp, lfm, lfc, lfp
                , clm, clc, clp
                , dxm, dxp, fdm, fdp, fdms, fdps
                , pm, pc, pp
                , edm, edc, edp
                , i, j, k, Comp::v());
      A.s[i][j][k] = tscn * am * aflagm;
      A.c[i][j][k]+= tscn * ac;
      A.n[i][j][k] = tscn * ap * aflagp;
  
      /* k-direction */
      onm=dom->ibody().on(i,j,k-1);
      onc=dom->ibody().on(i,j,k  );
      onp=dom->ibody().on(i,j,k+1);
      ofm=dom->ibody().off(i,j,k-1);
      ofc=dom->ibody().off(i,j,k  );
      ofp=dom->ibody().off(i,j,k+1);
      lsm=solid()->lambda(i,j,k-1);
      lsc=solid()->lambda(i,j,k  );
      lsp=solid()->lambda(i,j,k+1);
      lfm=fluid()->lambda(i,j,k-1);
      lfc=fluid()->lambda(i,j,k  );
      lfp=fluid()->lambda(i,j,k+1);
      clm=(*clr)[i][j][k-1];
      clc=(*clr)[i][j][k  ];
      clp=(*clr)[i][j][k+1];
      dxm=phi.dzb(k);
      dxp=phi.dzt(k);
      fdm=dom->ibody().fdzb(i,j,k);
      fdp=dom->ibody().fdzt(i,j,k);
      fdms=dom->ibody().fdzt(i,j,k-1);
      fdps=dom->ibody().fdzb(i,j,k+1);
      pm=phi[i][j][k-1];
      pc=phi[i][j][k  ];
      pp=phi[i][j][k+1];
      edm=edc=edp=0.0;
      if(diff_eddy){
        if (onm) edm = (*diff_eddy)[i][j][k-1];
        if (onc) edc = (*diff_eddy)[i][j][k  ];
        if (onp) edp = (*diff_eddy)[i][j][k+1];
      }
  
      diff_matrix(am, ac, ap
                , tm, tc, tp
                , aflagm, aflagp
                , vol, dSz(i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp, lfm, lfc, lfp
                , clm, clc, clp
                , dxm, dxp, fdm, fdp, fdms, fdps
                , pm, pc, pp
                , edm, edc, edp
                , i, j, k, Comp::w());
      A.b[i][j][k] = tscn * am * aflagm;
      A.c[i][j][k]+= tscn * ac;
      A.t[i][j][k] = tscn * ap * aflagp;
    }
  } /* conduction through solid */
 
  /*---------------------+ 
  |  this is needed too  |
  +---------------------*/
  if( !solid() ) 
    if(dom->ibody().nccells() > 0) {
      for_ijk(i,j,k)
        if( dom->ibody().off(i,j,k) ) {
          A.c[i][j][k]  = 1.0e-12;
          A.w[i][j][k]  = 0.0;
          A.e[i][j][k]  = 0.0;
          A.s[i][j][k]  = 0.0;
          A.n[i][j][k]  = 0.0;
          A.b[i][j][k]  = 0.0;
          A.t[i][j][k]  = 0.0;
          A.ci[i][j][k] = 1.0;
          fnew[i][j][k] = 0.0;
        }
    } /* is there an immersed body */

  A.c.exchange();
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();


#if 0
   std::cout<<"system_diffusive: end.\n";
#endif

}
