#include "enthalpytif.h"
#define VERSION_STABLE
using namespace std;

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void EnthalpyTIF::convection() {
  convection(&cnew);
}

/***************************************************************************//**
*  \brief Computes convection term using the last available velocities. 
*
*  \param conv - array into which convection term will be stored,
*
*  If called from new_time_step(), it will create \f$ \{C\}^{N-1} \f$; 
*  if called from inner iteration loop (from SIMPLE algorithm) it will compute
*  \f$ \{C\}^{N} \f$. 
*******************************************************************************/
void EnthalpyTIF::convection(Scalar * conv) {

  phi.exchange();

  real phim, phip;
  real umf, upf, vmf, vpf, wmf, wpf;

  for_aijk(i,j,k)
    (*conv)[i][j][k] = 0.0;

  /* find where are dirichlet or inlet boundary conditions */
  bool imin, imax, jmin, jmax, kmin, kmax; 
  imin = imax = jmin = jmax = kmin = kmax = false;
  for( int b=0; b<phi.bc().count(); b++ ) {
    if( phi.bc().type(b)==BndType::dirichlet() || 
        phi.bc().type(b)==BndType::inlet()     ||
        phi.bc().type(b)==BndType::outlet()    ||
        phi.bc().type(b)==BndType::neumann()   ||
        phi.bc().type(b)==BndType::insert() ) {
      if(phi.bc().direction(b)==Dir::imin() && dom->coord(Comp::i())==0)
        imin=true;
      if(phi.bc().direction(b)==Dir::imax() &&
         dom->coord(Comp::i())==dom->dim(Comp::i())-1) 
        imax=true;
      if(phi.bc().direction(b)==Dir::jmin() && dom->coord(Comp::j())==0)
        jmin=true;
      if(phi.bc().direction(b)==Dir::jmax() &&
         dom->coord(Comp::j())==dom->dim(Comp::j())-1) 
        jmax=true;
      if(phi.bc().direction(b)==Dir::kmin() && dom->coord(Comp::k())==0)
        kmin=true;
      if(phi.bc().direction(b)==Dir::kmax() &&
         dom->coord(Comp::k())==dom->dim(Comp::k())-1) 
        kmax=true;
    }
  }

  //boil::plot->plot(*clr, phi, iflag,"clr-phi-iflag", time->current_step());

  /* set flag */
  setflag();

  /*-----------------------+
  |  finite volume method  |
  +-----------------------*/
  for_ijk(i,j,k) {
    { /////////
      //  u  //
      /////////
    umf = (*u)[Comp::u()][i]  [j][k];  // u @ imin
    upf = (*u)[Comp::u()][i+1][j][k];  // u @ imax
    
    real a_w = dSx(i,j,k);
    real a_e = dSx(i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_w *= dom->ibody().fSw(i,j,k);
      a_e *= dom->ibody().fSe(i,j,k);
    }

    phim = a_w*lim.limit(-umf, phi[i+1][j][k],phi[i][j][k],phi[i-1][j][k]);
    phip = a_e*lim.limit(+upf, phi[i-1][j][k],phi[i][j][k],phi[i+1][j][k]);
  
    if(i==si() && imin) phim = a_w*phi[i-1][j][k];
    if(i==ei() && imax) phip = a_e*phi[i+1][j][k];

    (*conv)[i]  [j][k] += (umf*phim - upf*phip);
    (*conv)[i+1][j][k] +=  upf*phip;
    (*conv)[i-1][j][k] -=  umf*phim;
    }
    { /////////
      //  v  //
      /////////
    vmf = (*u)[Comp::v()][i][j]  [k];  // v @ jmin
    vpf = (*u)[Comp::v()][i][j+1][k];  // v @ jmax

    real a_s = dSy(i,j,k);
    real a_n = dSy(i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_s *= dom->ibody().fSs(i,j,k);
      a_n *= dom->ibody().fSn(i,j,k);
    }

    phim = a_s*lim.limit(-vmf, phi[i][j+1][k],phi[i][j][k],phi[i][j-1][k]);
    phip = a_n*lim.limit(+vpf, phi[i][j-1][k],phi[i][j][k],phi[i][j+1][k]);

    if(j==sj() && jmin) phim = a_s*phi[i][j-1][k];
    if(j==ej() && jmax) phip = a_n*phi[i][j+1][k];

    (*conv)[i][j]  [k] += (vmf*phim - vpf*phip);
    (*conv)[i][j+1][k] +=  vpf*phip;
    (*conv)[i][j-1][k] -=  vmf*phim;
    } 
    { /////////
      //  w  //
      /////////
    wmf = (*u)[Comp::w()][i][j][k];    // w @ kmin
    wpf = (*u)[Comp::w()][i][j][k+1];  // w @ kmax

    real a_b = dSz(i,j,k);
    real a_t = dSz(i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_b *= dom->ibody().fSb(i,j,k);
      a_t *= dom->ibody().fSt(i,j,k);
    }

    phim = a_b*lim.limit(-wmf, phi[i][j][k+1], phi[i][j][k], phi[i][j][k-1]);
    phip = a_t*lim.limit(+wpf, phi[i][j][k-1], phi[i][j][k], phi[i][j][k+1]);

    if(k==sk() && kmin) phim = a_b*phi[i][j][k-1];
    if(k==ek() && kmax) phip = a_t*phi[i][j][k+1];

    (*conv)[i][j][k]   += (wmf*phim - wpf*phip);
    (*conv)[i][j][k+1] +=  wpf*phip;
    (*conv)[i][j][k-1] -=  wmf*phim;
    }
  }

  // k-direction
  for_ij(i,j) buff[i][j][ek()] = (*conv)[i][j][ek()+1]; 
  for_ij(i,j) buff[i][j][sk()] = (*conv)[i][j][sk()-1]; 
  buff.exchange(2);
  for_ij(i,j) (*conv)[i][j][ek()] += buff[i][j][ek()+1]; 
  for_ij(i,j) (*conv)[i][j][sk()] += buff[i][j][sk()-1]; 

  // j-direction
  for_ik(i,k) buff[i][ej()][k] = (*conv)[i][ej()+1][k]; 
  for_ik(i,k) buff[i][sj()][k] = (*conv)[i][sj()-1][k]; 
  buff.exchange(1);
  for_ik(i,k)  (*conv)[i][ej()][k] += buff[i][ej()+1][k]; 
  for_ik(i,k)  (*conv)[i][sj()][k] += buff[i][sj()-1][k]; 

  // i-direction
  for_jk(j,k) buff[ei()][j][k] = (*conv)[ei()+1][j][k]; 
  for_jk(j,k) buff[si()][j][k] = (*conv)[si()-1][j][k]; 
  buff.exchange(0);
  for_jk(j,k) (*conv)[ei()][j][k] += buff[ei()+1][j][k];
  for_jk(j,k) (*conv)[si()][j][k] += buff[si()-1][j][k];

  for_ijk(i,j,k) {
    real divu = - dSx(i,j,k)*(*u)[Comp::u()][i]  [j]  [k]
                + dSx(i,j,k)*(*u)[Comp::u()][i+1][j]  [k]
                - dSy(i,j,k)*(*u)[Comp::v()][i]  [j]  [k]
                + dSy(i,j,k)*(*u)[Comp::v()][i]  [j+1][k]
                - dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k]
                + dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k+1];
    (*conv)[i][j][k] += phi[i][j][k] * divu;
    //(*conv)[i][j][k] /= dV(i,j,k);
  }

#if 1
  /*---------------------------+
  |  finite difference method  |
  +---------------------------*/
  for_ijk(i,j,k) {

    if(iflag[i][j][k]==0){

      real umf, upf, vmf, vpf, wmf, wpf, uc, vc, wc;
      real dtdxm, dtdxp, dtdym, dtdyp, dtdzm, dtdzp;
      real udtdx, vdtdy, wdtdz;

      real clrc = clrold[i][j][k];

      /* set velocity */
      int iumf,iupf,ivmf,ivpf,iwmf,iwpf;
      iumf=iupf=ivmf=ivpf=iwmf=iwpf=0;

      // u
      umf = (*u)[Comp::u()][i]  [j][k];
      upf = (*u)[Comp::u()][i+1][j][k];

      // v
      vmf = (*u)[Comp::v()][i][j]  [k];
      vpf = (*u)[Comp::v()][i][j+1][k];

      // w
      wmf = (*u)[Comp::w()][i][j][k];
      wpf = (*u)[Comp::w()][i][j][k+1];

      // dtdxm
      real clrw = clrold[i-1][j][k];
      //if((clrc-clrsurf)*(clrw-clrsurf)<0.0){
      if(Interface_old(-1,Comp::i(),i,j,k)) {
        real dxm, ts;
        if(!fs) {
          dxm = max(epsl,(clrsurf-clrc)/(clrw-clrc));
          ts = Tint_old(-1,Comp::i(),dxm,i,j,k);
          dxm = dxm * dxw(i);
        } else {
          dxm = max(epsl*dxw(i),distance_x(i,j,k,-1,ts,true));
        }
        dtdxm = (phi[i][j][k]-ts)/dxm;
#ifndef VERSION_STABLE
        umf = upf;
#endif
      } else {
        if(dom->ibody().off(i-1,j,k)) {
          dtdxm = gradt_ib(-1,Comp::i(),i,j,k);
        } else {
          dtdxm = (phi[i  ][j][k]-phi[i-1][j][k])/dxw(i);
        }
      }

      // dtdxp
      real clre = clrold[i+1][j][k];
      //if((clrc-clrsurf)*(clre-clrsurf)<0.0){
      if(Interface_old(+1,Comp::i(),i,j,k)) {
        real dxp, ts;
        if(!fs) {
          dxp = max(epsl,(clrsurf-clrc)/(clre-clrc));
          ts = Tint_old(+1,Comp::i(),dxp,i,j,k);
          dxp = dxp * dxe(i);
        } else {
          dxp = max(epsl*dxe(i),distance_x(i,j,k,+1,ts,true));
        }
        dtdxp = (ts-phi[i][j][k])/dxp;
#ifndef VERSION_STABLE
        upf = umf;
#endif
      } else {
        if(dom->ibody().off(i+1,j,k)) {
          dtdxp = gradt_ib(+1,Comp::i(),i,j,k);
        } else {
          dtdxp = (phi[i+1][j][k]-phi[i  ][j][k])/dxe(i);
        }
      }

      // dtdym
      real clrs = clrold[i][j-1][k];
      //if((clrc-clrsurf)*(clrs-clrsurf)<0.0){
      if(Interface_old(-1,Comp::j(),i,j,k)) {
        real dym, ts;
        if(!fs) {
          dym = max(epsl,(clrsurf-clrc)/(clrs-clrc));
          ts = Tint_old(-1,Comp::j(),dym,i,j,k);
          dym = dym * dys(j);
        } else {
          dym = max(epsl*dys(j),distance_y(i,j,k,-1,ts,true));
        }
        dtdym = (phi[i][j][k]-ts)/dym;
#ifndef VERSION_STABLE
        vmf = vpf;
#endif
      } else {
        if(dom->ibody().off(i,j-1,k)) {
          dtdym = gradt_ib(-1,Comp::j(),i,j,k);
        } else {
          dtdym = (phi[i][j  ][k]-phi[i][j-1][k])/dys(j);
        }
      }

      // dtdyp
      real clrn = clrold[i][j+1][k];
      //if((clrc-clrsurf)*(clrn-clrsurf)<0.0){
      if(Interface_old(+1,Comp::j(),i,j,k)) {
        real dyp, ts;
        if(!fs) {
          dyp = max(epsl,(clrsurf-clrc)/(clrn-clrc));
          ts = Tint_old(+1,Comp::j(),dyp,i,j,k);
          dyp = dyp * dyn(j);
        } else {
          dyp = max(epsl*dyn(j),distance_y(i,j,k,+1,ts,true));
        }
        dtdyp = (ts-phi[i][j][k])/dyp;
#ifndef VERSION_STABLE
        vpf = vmf;
#endif
      } else {
        if(dom->ibody().off(i,j+1,k)) {
          dtdyp = gradt_ib(+1,Comp::j(),i,j,k);
        } else {
          dtdyp = (phi[i][j+1][k]-phi[i][j  ][k])/dyn(j);
        }
      }

      // dtdzm
      real clrb = clrold[i][j][k-1];
      //if((clrc-clrsurf)*(clrb-clrsurf)<0.0){
      if(Interface_old(-1,Comp::k(),i,j,k)) {
        real dzm, ts;
        if(!fs) {
          dzm = max(epsl,(clrsurf-clrc)/(clrb-clrc));
          ts = Tint_old(-1,Comp::k(),dzm,i,j,k);
          dzm = dzm * dzb(k);
        } else {
          dzm = max(epsl*dzb(k),distance_z(i,j,k,-1,ts,true));
        }
        dtdzm = (phi[i][j][k]-ts)/dzm;
#ifndef VERSION_STABLE
        wmf = wpf;
#endif
      } else {
        if(dom->ibody().off(i,j,k-1)) {
          dtdzm = gradt_ib(-1,Comp::k(),i,j,k);
        } else {
          dtdzm = (phi[i][j][k  ]-phi[i][j][k-1])/dzb(k);
        }
      }

      // dtdzp
      real clrt = clrold[i][j][k+1];
      //if((clrc-clrsurf)*(clrt-clrsurf)<0.0){
      if(Interface_old(+1,Comp::k(),i,j,k)) {
        real dzp, ts;
        if(!fs) {
          dzp = max(epsl,(clrsurf-clrc)/(clrt-clrc));
          ts = Tint_old(+1,Comp::k(),dzp,i,j,k);
          dzp = dzp * dzt(k);
        } else {
          dzp = max(epsl*dzt(k),distance_z(i,j,k,+1,ts,true));
        }
        dtdzp = (ts-phi[i][j][k])/dzp;
#ifndef VERSION_STABLE
        wpf = wmf;
#endif
      } else {
        if(dom->ibody().off(i,j,k+1)) {
          dtdzp = gradt_ib(+1,Comp::k(),i,j,k);
        } else {
          dtdzp = (phi[i][j][k+1]-phi[i][j][k  ])/dzt(k);
        }
      }

      uc  = 0.5*(umf+upf);
      vc  = 0.5*(vmf+vpf);
      wc  = 0.5*(wmf+wpf);
#if 1
      udtdx = 0.5*(uc+fabs(uc))*dtdxm
            + 0.5*(uc-fabs(uc))*dtdxp;
      vdtdy = 0.5*(vc+fabs(vc))*dtdym
            + 0.5*(vc-fabs(vc))*dtdyp;
      wdtdz = 0.5*(wc+fabs(wc))*dtdzm
            + 0.5*(wc-fabs(wc))*dtdzp;
#else
      real u_sign=copysign(1.0,uc);
      real v_sign=copysign(1.0,vc);
      real w_sign=copysign(1.0,wc);
      udtdx = 0.5 * (u_sign+fabs(u_sign)) * dtdxm * umf  // uc > 0
            - 0.5 * (u_sign-fabs(u_sign)) * dtdxp * upf; // uc < 0
      vdtdy = 0.5 * (v_sign+fabs(v_sign)) * dtdym * vmf  // vc > 0
            - 0.5 * (v_sign-fabs(v_sign)) * dtdyp * vpf; // vc < 0
      wdtdz = 0.5 * (w_sign+fabs(w_sign)) * dtdzm * wmf  // wc > 0
            - 0.5 * (w_sign-fabs(w_sign)) * dtdzp * wpf; // wc < 0

#endif
      (*conv)[i][j][k] = - (udtdx + vdtdy + wdtdz) * dV(i,j,k);

    }
  }
#endif

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k){
      if(dom->ibody().off(i,j,k)){
        (*conv)[i][j][k] = 0.0;
      }
    }
  }

  for_ijk(i,j,k) {
    real r,c;
    if(clrold[i][j][k]>=clrsurf){
      c = cpl;
    } else {
      c = cpv;
    }
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }

}

