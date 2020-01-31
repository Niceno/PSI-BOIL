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
  if( accelerated_no_solid ) {

    /* coefficients in i direction (w and e) */
    for_ijk(i,j,k) {
      real lc, cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
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
      real xm,xp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if(!Interface(-1,Comp::i(),i,j,k)){
        xm=phi.dxw(i);
      } else {
        real ts;
        xm = std::max(epsl*phi.dxw(i),distance_x(i,j,k,-1,ts));
        aflagm=0.0;
      }
      if(!Interface(+1,Comp::i(),i,j,k)){
        xp=phi.dxe(i);
      } else {
        real ts;
        xp = std::max(epsl*phi.dxe(i),distance_x(i,j,k,+1,ts));
        aflagp=0.0;
      }
      real cxm = coef_x_m(xm,xp,phi.xc(i));
      real cxp = coef_x_p(xm,xp,phi.xc(i));

      A.w[i][j][k] =  tscn * lc * vol * cxm * aflagm;
      A.c[i][j][k] += tscn * lc * vol * (cxm+cxp);
      A.e[i][j][k] =  tscn * lc * vol * cxp * aflagp;
    }

    /* coefficients in j direction (s and n) */
    for_ijk(i,j,k) {
      real lc, cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
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
      real ym,yp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if(!Interface(-1,Comp::j(),i,j,k)){
        ym=phi.dys(j);
      } else {
        real ts;
        ym = std::max(epsl*phi.dys(j),distance_y(i,j,k,-1,ts));
        aflagm=0.0;
      }
      if(!Interface(+1,Comp::j(),i,j,k)){
        yp=phi.dyn(j);
      } else {
        real ts;
        yp = std::max(epsl*phi.dyn(j),distance_y(i,j,k,+1,ts));
        aflagp=0.0;
      }
      real cym = coef_y_m(ym,yp,phi.yc(j));
      real cyp = coef_y_p(ym,yp,phi.yc(j));

      A.s[i][j][k] =  tscn * lc * vol * cym * aflagm;
      A.c[i][j][k] += tscn * lc * vol * (cym+cyp);
      A.n[i][j][k] =  tscn * lc * vol * cyp * aflagp;
    }

    /* coefficients in k direction (b and t) */
    for_ijk(i,j,k) {
      real lc, cp_mass;
      if((*clr)[i][j][k]>=clrsurf){
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
      real zm,zp,aflagm,aflagp;
      aflagm=aflagp=1.0;
      if(!Interface(-1,Comp::k(),i,j,k)){
        zm=phi.dzb(k);
      } else {
        real ts;
        zm = std::max(epsl*phi.dzb(k),distance_z(i,j,k,-1,ts));
        aflagm=0.0;
      }
      if(!Interface(+1,Comp::k(),i,j,k)){
        zp=phi.dzt(k);
      } else {
        real ts;
        zp = std::max(epsl*phi.dzt(k),distance_z(i,j,k,+1,ts));
        aflagp=0.0;
      }
      real czm = coef_z_m(zm,zp,phi.zc(k));
      real czp = coef_z_p(zm,zp,phi.zc(k));

      A.b[i][j][k] =  tscn * lc * vol * czm * aflagm;
      A.c[i][j][k] += tscn * lc * vol * (czm+czp);
      A.t[i][j][k] =  tscn * lc * vol * czp * aflagp;
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

    for_ijk(i,j,k) {

      const real vol = phi.dV(i,j,k);
      real lc, xm, xp;

      bool onm, onc, onp, ofm, ofc, ofp; // on & off
      real lsm, lsc, lsp; // lambda
      int clm, clc, clp; // topology flag
      real dxm, dxp, fdm, fdp, fdms, fdps;
      real pm, pc, pp;    // temperature-input
      real edm, edc, edp; // eddy viscosity
      real am, ac, ap;
      real tm, tc, tp;    // temperature-output
      real aflagm, aflagp;
      real pos0;
      coef_gen coef_m, coef_p;

      /* i-direction */
      onm=dom->ibody().on(i-1,j,k);
      onc=dom->ibody().on(i  ,j,k);
      onp=dom->ibody().on(i+1,j,k);
      ofm=dom->ibody().off(i-1,j,k);
      ofc=dom->ibody().off(i  ,j,k);
      ofp=dom->ibody().off(i+1,j,k);
      lsm=safe_solid->lambda(i-1,j,k);
      lsc=safe_solid->lambda(i  ,j,k);
      lsp=safe_solid->lambda(i+1,j,k);
      clm=iflag[i-1][j][k];
      clc=iflag[i  ][j][k];
      clp=iflag[i+1][j][k];
 
      dxm=phi.dxw(i);
      dxp=phi.dxe(i);
      pos0=phi.xc(i);
      coef_m = &EnthalpyFD::coef_x_m;
      coef_p = &EnthalpyFD::coef_x_p;
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
                , pos0, coef_m, coef_p
                , vol, dSx(Sign::neg(),i,j,k), dSx(Sign::pos(),i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp
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
      lsm=safe_solid->lambda(i,j-1,k);
      lsc=safe_solid->lambda(i,j  ,k);
      lsp=safe_solid->lambda(i,j+1,k);
      clm=iflag[i][j-1][k];
      clc=iflag[i][j  ][k];
      clp=iflag[i][j+1][k];
      dxm=phi.dys(j);
      dxp=phi.dyn(j);
      pos0=phi.yc(j);
      coef_m = &EnthalpyFD::coef_y_m;
      coef_p = &EnthalpyFD::coef_y_p;
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
                , pos0, coef_m, coef_p
                , vol, dSy(Sign::neg(),i,j,k), dSy(Sign::pos(),i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp
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
      lsm=safe_solid->lambda(i,j,k-1);
      lsc=safe_solid->lambda(i,j,k  );
      lsp=safe_solid->lambda(i,j,k+1);
      clm=iflag[i][j][k-1];
      clc=iflag[i][j][k  ];
      clp=iflag[i][j][k+1];
      dxm=phi.dzb(k);
      dxp=phi.dzt(k);
      pos0=phi.zc(k);
      coef_m = &EnthalpyFD::coef_z_m;
      coef_p = &EnthalpyFD::coef_z_p;
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
                , pos0, coef_m, coef_p
                , vol, dSz(Sign::neg(),i,j,k), dSz(Sign::pos(),i,j,k)
                , onm, onc, onp, ofm, ofc, ofp
                , lsm, lsc, lsp
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
