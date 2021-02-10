#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$
*         or adds diffusion to right hand side.
*******************************************************************************/
void EnthalpyFD::evaluate_diffusion(const Old old, const Scalar * diff_eddy) {

  /* initialize: get time stepping coefficient */
  real tscn = diff_ts.N(); /* 1.0 for fully implicit */
  real tscm =  diff_ts.Nm1(); /* 0.5 for c.n. 0.0 for fully implicit */

  ResistEval re;
  if(old==Old::no) {
    assert( tscn > 0.0 );
    re = ResistEval::no;
    ftif = 0.;
  } else {
    assert( tscm > 0.0 );
    re = ResistEval::yes;
  }

  /*------------------------------------+
  |  no conduction through solid parts  |
  +------------------------------------*/
  if( accelerated_no_solid ) {

    std::array<ConnectType,3> ctype = { ConnectType::fluid, 
                                        ConnectType::fluid,
                                        ConnectType::fluid };

    std::vector<StencilPoint> stencil;

    /* coefficients in i direction (w and e) */
    for_ijk(i,j,k) {
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      real xm,xp,pm,pp;
      stencil.clear();
      if(!cht.interface(Sign::neg(),Comp::i(),i,j,k)){
        xm=phi.dxw(i);
        pm=phi[i-1][j][k];
        ctype[0] = ConnectType::fluid;
      } else {
        xm = cht.distance_int_x(Sign::neg(),i,j,k,pm,re,old);
        ctype[0] = ConnectType::interface;
      }
      if(!cht.interface(Sign::pos(),Comp::i(),i,j,k)){
        xp=phi.dxe(i);
        pp=phi[i+1][j][k];
        ctype[2] = ConnectType::fluid;
      } else {
        xp = cht.distance_int_x(Sign::pos(),i,j,k,pp,re,old);
        ctype[2] = ConnectType::interface;
      }
      real cxm = coef_x_m(xm,xp,phi.xc(i)) * lc * vol;
      real cxp = coef_x_p(xm,xp,phi.xc(i)) * lc * vol;

      if(old==Old::no) {
        stencil.push_back(StencilPoint(0,pm,-xm));
        stencil.push_back(StencilPoint(1,phi[i][j][k],0.));
        stencil.push_back(StencilPoint(2,pp, xp));
        diffmatrix_kernel(ctype,tscn*cxm,tscn*cxp,stencil,
                          A.w[i][j][k],A.c[i][j][k],A.e[i][j][k],ftif[i][j][k]);
      } else {
        const real Aw = tscm * cxm;
        const real Ac = tscm * (cxm+cxp);
        const real Ae = tscm * cxp;
        fold[i][j][k] += Aw*pm - Ac*phi[i][j][k] + Ae*pp;
      }
    }

    /* coefficients in j direction (s and n) */
    for_ijk(i,j,k) {
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      real ym,yp,pm,pp;
      stencil.clear();
      if(!cht.interface(Sign::neg(),Comp::j(),i,j,k)){
        ym=phi.dys(j);
        pm=phi[i][j-1][k];
        ctype[0] = ConnectType::fluid;
      } else {
        ym = cht.distance_int_y(Sign::neg(),i,j,k,pm,re,old);
        ctype[0] = ConnectType::interface;
      }
      if(!cht.interface(Sign::pos(),Comp::j(),i,j,k)){
        yp=phi.dyn(j);
        pp=phi[i][j+1][k];
        ctype[2] = ConnectType::fluid;
      } else {
        yp = cht.distance_int_y(Sign::pos(),i,j,k,pp,re,old);
        ctype[2] = ConnectType::interface;
      }
      real cym = coef_y_m(ym,yp,phi.yc(j)) * lc * vol;
      real cyp = coef_y_p(ym,yp,phi.yc(j)) * lc * vol;

      if(old==Old::no) {
        stencil.push_back(StencilPoint(0,pm,-ym));
        stencil.push_back(StencilPoint(1,phi[i][j][k],0.));
        stencil.push_back(StencilPoint(2,pp, yp));
        diffmatrix_kernel(ctype,tscn*cym,tscn*cyp,stencil,
                          A.s[i][j][k],A.c[i][j][k],A.n[i][j][k],ftif[i][j][k]);
      } else {
        const real As = tscm * cym;
        const real Ac = tscm * (cym+cyp);
        const real An = tscm * cyp;
        fold[i][j][k] += As*pm - Ac*phi[i][j][k] + An*pp;
      }
    }

    /* coefficients in k direction (b and t) */
    for_ijk(i,j,k) {
      real lc = cht.lambda(i,j,k,diff_eddy);
      const real vol = phi.dV(i,j,k);
      real zm,zp,pm,pp;
      stencil.clear();
      if(!cht.interface(Sign::neg(),Comp::k(),i,j,k)){
        zm=phi.dzb(k);
        pm=phi[i][j][k-1];
        ctype[0] = ConnectType::fluid;
      } else {
        zm = cht.distance_int_z(Sign::neg(),i,j,k,pm,re,old);
        ctype[0] = ConnectType::interface;
      }
      if(!cht.interface(Sign::pos(),Comp::k(),i,j,k)){
        zp=phi.dzt(k);
        pp=phi[i][j][k+1];
        ctype[2] = ConnectType::fluid;
      } else {
        zp = cht.distance_int_z(Sign::pos(),i,j,k,pp,re,old);
        ctype[2] = ConnectType::interface;
      }
      real czm = coef_z_m(zm,zp,phi.zc(k)) * lc * vol;
      real czp = coef_z_p(zm,zp,phi.zc(k)) * lc * vol;

      if(old==Old::no) {
        stencil.push_back(StencilPoint(0,pm,-zm));
        stencil.push_back(StencilPoint(1,phi[i][j][k],0.));
        stencil.push_back(StencilPoint(2,pp, zp));
        diffmatrix_kernel(ctype,tscn*czm,tscn*czp,stencil,
                          A.b[i][j][k],A.c[i][j][k],A.t[i][j][k],ftif[i][j][k]);
      } else {
        const real Ab = tscm * czm;
        const real Ac = tscm * (czm+czp);
        const real At = tscm * czp;
        fold[i][j][k] += Ab*pm - Ac*phi[i][j][k] + At*pp;
      }
    }
    /*-------------------------------+
    |  a "touch" from immersed body  |
    +-------------------------------*/
    if(dom->ibody().nccells() > 0) {
      boil::oout<<"EFD: Conduction without solid. "
                <<"Underdevelopment. Exiting."<<boil::endl;
      exit(0); 
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

    for_m(m) {
      int ii(0),jj(0),kk(0);
      if       (m==Comp::i()) {
        ii=1;
      } else if(m==Comp::j()) {
        jj=1;
      } else {
        kk=1;
      }

      for_ijk(i,j,k) {

        const real vol = phi.dV(i,j,k);
        bool onm, onc, onp, ofm, ofc, ofp; /* on & off */
        real lsm, lsc, lsp; /* lambdas */
        real lvm, lvc, lvp; /* lambdav */
        real llm, llc, llp; /* lambdal */
        int clm, clc, clp;  /* topology flag */
        real dxm, dxp, fdm, fdp, fdms, fdps;
        real am, ac, ap;
        real tm, tc, tp;
        bool aflagm, aflagp;
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

        tm=phi[i-ii][j-jj][k-kk];
        tc=phi[i   ][j   ][k   ];
        tp=phi[i+ii][j+jj][k+kk];

        if       (m==Comp::i()) {
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
        } else if(m==Comp::j()) {
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
                  , i, j, k, m);

        if(old==Old::no) {
          if       (m==Comp::i()) {
            A.w[i][j][k] = aflagm ? tscn * am : 0.0;
            A.c[i][j][k]+= tscn * ac;
            A.e[i][j][k] = aflagp ? tscn * ap : 0.0;
          } else if(m==Comp::j()) {
            A.s[i][j][k] = aflagm ? tscn * am : 0.0;
            A.c[i][j][k]+= tscn * ac;
            A.n[i][j][k] = aflagp ? tscn * ap : 0.0;
          } else {
            A.b[i][j][k] = aflagm ? tscn * am : 0.0;
            A.c[i][j][k]+= tscn * ac;
            A.t[i][j][k] = aflagp ? tscn * ap : 0.0;
          }

          /* rhs */
          ftif[i][j][k] += tscn * sourceterm;
          ftif[i][j][k] += aflagm ? 0.0 : tscn * am * tm;
          ftif[i][j][k] += aflagp ? 0.0 : tscn * ap * tp;
        } else {
          fold[i][j][k] -= tscm * ac*tc;
          fold[i][j][k] += tscm * am*tm;
          fold[i][j][k] += tscm * ap*tp;
          fold[i][j][k] += tscm * sourceterm;
        }
     
      } /* ijk */
    } /* m */

  } /* conduction through solid */
 
  if(old==Old::no) {
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
  }

  return;
}
