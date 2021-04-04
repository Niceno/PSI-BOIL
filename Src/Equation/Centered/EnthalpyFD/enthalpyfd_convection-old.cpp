#include "enthalpyfd.h"
#include "def.h"

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void EnthalpyFD::convection() {
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
void EnthalpyFD::convection(Scalar * conv) {

#ifdef CNEW
  Old old = Old::no;
#else
  Old old = Old::yes;
#endif

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

  //boil::plot->plot(*clr, phi, iflagold,"clr-phi-iflag", time->current_step());

  for_ijk(i,j,k) {
    { /////////
      //  u  //
      /////////
    if(cht.above_interface(i,j,k,old)) {
      umf = (*uliq)[Comp::u()][i]  [j][k];  // u @ imin
      upf = (*uliq)[Comp::u()][i+1][j][k];  // u @ imax
    } else {
      umf = (*ugas)[Comp::u()][i]  [j][k];  // u @ imin
      upf = (*ugas)[Comp::u()][i+1][j][k];  // u @ imax
    }
    
    if(umf>0.) {
      phim = lim.limit(+umf, phi[i-2][j][k],phi[i-1][j][k],phi[i][j][k]);
    } else {
      phim = lim.limit(-umf, phi[i+1][j][k],phi[i][j][k],phi[i-1][j][k]);
    }
    if(upf>0.) {
      phip = lim.limit(+upf, phi[i-1][j][k],phi[i][j][k],phi[i+1][j][k]);
    } else {
      phip = lim.limit(-upf, phi[i+2][j][k],phi[i+1][j][k],phi[i][j][k]);
    }
  
    if(i==si() && imin) phim = phi[i-1][j][k];
    if(i==ei() && imax) phip = phi[i+1][j][k];

    flux_liq[Comp::u()][i][j][k] = umf*phim;
    flux_liq[Comp::u()][i+1][j][k] = upf*phip;

    }
    { /////////
      //  v  //
      /////////
    if(cht.above_interface(i,j,k,old)) {
      vmf = (*uliq)[Comp::v()][i][j]  [k];  // v @ jmin
      vpf = (*uliq)[Comp::v()][i][j+1][k];  // v @ jmax
    } else {
      vmf = (*ugas)[Comp::v()][i][j]  [k];  // v @ jmin
      vpf = (*ugas)[Comp::v()][i][j+1][k];  // v @ jmax
    }

    if(vmf>0.) {
      phim = lim.limit(+vmf, phi[i][j-2][k],phi[i][j-1][k],phi[i][j][k]);
    } else {
      phim = lim.limit(-vmf, phi[i][j+1][k],phi[i][j][k],phi[i][j-1][k]);
    }
    if(vpf>0.) {
      phip = lim.limit(+vpf, phi[i][j-1][k],phi[i][j][k],phi[i][j+1][k]);
    } else {
      phip = lim.limit(-vpf, phi[i][j+2][k],phi[i][j+1][k],phi[i][j][k]);
    }

    if(j==sj() && jmin) phim = phi[i][j-1][k];
    if(j==ej() && jmax) phip = phi[i][j+1][k];

    flux_liq[Comp::v()][i][j][k] = vmf*phim;
    flux_liq[Comp::v()][i][j+1][k] = vpf*phip;

    } 
    { /////////
      //  w  //
      /////////
    if(cht.above_interface(i,j,k,old)) {
      wmf = (*uliq)[Comp::w()][i][j][k];    // w @ kmin
      wpf = (*uliq)[Comp::w()][i][j][k+1];  // w @ kmax
    } else {
      wmf = (*ugas)[Comp::w()][i][j][k];    // w @ kmin
      wpf = (*ugas)[Comp::w()][i][j][k+1];  // w @ kmax
    }

    if(wmf>0.) {
      phim = lim.limit(+wmf, phi[i][j][k-2],phi[i][j][k-1],phi[i][j][k]);
    } else {
      phim = lim.limit(-wmf, phi[i][j][k+1], phi[i][j][k], phi[i][j][k-1]);
    }
    if(wpf>0.) {
      phip = lim.limit(+wpf, phi[i][j][k-1], phi[i][j][k], phi[i][j][k+1]);
    } else {
      phip = lim.limit(-wpf, phi[i][j][k+2],phi[i][j][k+1],phi[i][j][k]);
    }

    if(k==sk() && kmin) phim = phi[i][j][k-1];
    if(k==ek() && kmax) phip = phi[i][j][k+1];

    flux_liq[Comp::w()][i][j][k] = wmf*phim;
    flux_liq[Comp::w()][i][j][k+1] = wpf*phip;

    }

    (*conv)[i][j][k] = -flux_liq.divergence(i,j,k)*dV(i,j,k);
  }

#if 1
  /*---------------------------+
  |  finite difference method  |
  +---------------------------*/
  for_ijk(i,j,k) {

  #ifdef USE_FDM_CONVECTION
    if(true) {
  #else
    #ifdef CNEW
    if(abs(iflag[i][j][k])<3){
    #else
    if(abs(iflagold[i][j][k])<3){
    #endif
  #endif
      real umf, upf, vmf, vpf, wmf, wpf, uc, vc, wc;
      real dtdxm, dtdxp, dtdym, dtdyp, dtdzm, dtdzp;
      real udtdx, vdtdy, wdtdz;

      /* set velocity */
      int iumf,iupf,ivmf,ivpf,iwmf,iwpf;
      iumf=iupf=ivmf=ivpf=iwmf=iwpf=0;

      if(cht.above_interface(i,j,k,old)) {
        // u
        umf = (*uliq)[Comp::u()][i]  [j][k];
        upf = (*uliq)[Comp::u()][i+1][j][k];

        // v
        vmf = (*uliq)[Comp::v()][i][j]  [k];
        vpf = (*uliq)[Comp::v()][i][j+1][k];

        // w
        wmf = (*uliq)[Comp::w()][i][j][k];
        wpf = (*uliq)[Comp::w()][i][j][k+1];
      } else {
        // u
        umf = (*ugas)[Comp::u()][i]  [j][k];
        upf = (*ugas)[Comp::u()][i+1][j][k];
  
        // v
        vmf = (*ugas)[Comp::v()][i][j]  [k];
        vpf = (*ugas)[Comp::v()][i][j+1][k];

        // w
        wmf = (*ugas)[Comp::w()][i][j][k];
        wpf = (*ugas)[Comp::w()][i][j][k+1];
      }

#if 0
  #ifdef CNEW
      dtdxm = dtdxp = cht.first_derivative(false,Comp::i(),i,j,k,
                                           ao_conv,false,Old::no);
      dtdym = dtdyp = cht.first_derivative(false,Comp::j(),i,j,k,
                                           ao_conv,false,Old::no);
      dtdzm = dtdzp = cht.first_derivative(false,Comp::k(),i,j,k,
                                           ao_conv,false,Old::no);
  #else
      dtdxm = dtdxp = cht.first_derivative(false,Comp::i(),i,j,k,
                                           ao_conv,false,Old::yes);
      dtdym = dtdyp = cht.first_derivative(false,Comp::j(),i,j,k,
                                           ao_conv,false,Old::yes);
      dtdzm = dtdzp = cht.first_derivative(false,Comp::k(),i,j,k,
                                           ao_conv,false,Old::yes);
  #endif
#else
      // dtdxm
      if(cht.interface(Sign::neg(),Comp::i(),i,j,k,old)) {
        real dxm, ts;
        dxm = cht.distance_int(Sign::neg(),Comp::i(),i,j,k,ts,ResistEval::yes,old);
        dtdxm = (phi[i][j][k]-ts)/dxm;
#ifndef VERSION_STABLE
        umf = upf;
#endif
      } else {
        if(dom->ibody().off(i-1,j,k)) {
          dtdxm = cht.gradt_ib(Sign::neg(),Comp::i(),i,j,k,old,phi);
        } else {
          dtdxm = (phi[i  ][j][k]-phi[i-1][j][k])/dxw(i);
        }
      }

      // dtdxp
      if(cht.interface(Sign::pos(),Comp::i(),i,j,k,old)) {
        real dxp, ts;
        dxp = cht.distance_int(Sign::pos(),Comp::i(),i,j,k,ts,ResistEval::yes,old);
        dtdxp = (ts-phi[i][j][k])/dxp;
#ifndef VERSION_STABLE
        upf = umf;
#endif
      } else {
        if(dom->ibody().off(i+1,j,k)) {
          dtdxp = cht.gradt_ib(Sign::pos(),Comp::i(),i,j,k,old,phi);
        } else {
          dtdxp = (phi[i+1][j][k]-phi[i  ][j][k])/dxe(i);
        }
      }

      // dtdym
      if(cht.interface(Sign::neg(),Comp::j(),i,j,k,old)) {
        real dym, ts;
        dym = cht.distance_int(Sign::neg(),Comp::j(),i,j,k,ts,ResistEval::yes,old);
        dtdym = (phi[i][j][k]-ts)/dym;
#ifndef VERSION_STABLE
        vmf = vpf;
#endif
      } else {
        if(dom->ibody().off(i,j-1,k)) {
          dtdym = cht.gradt_ib(Sign::neg(),Comp::j(),i,j,k,old,phi);
        } else {
          dtdym = (phi[i][j  ][k]-phi[i][j-1][k])/dys(j);
        }
      }

      // dtdyp
      if(cht.interface(Sign::pos(),Comp::j(),i,j,k,old)) {
        real dyp, ts;
        dyp = cht.distance_int(Sign::pos(),Comp::j(),i,j,k,ts,ResistEval::yes,old);
        dtdyp = (ts-phi[i][j][k])/dyp;
#ifndef VERSION_STABLE
        vpf = vmf;
#endif
      } else {
        if(dom->ibody().off(i,j+1,k)) {
          dtdyp = cht.gradt_ib(Sign::pos(),Comp::j(),i,j,k,old,phi);
        } else {
          dtdyp = (phi[i][j+1][k]-phi[i][j  ][k])/dyn(j);
        }
      }

      // dtdzm
      if(cht.interface(Sign::neg(),Comp::k(),i,j,k,old)) {
        real dzm, ts;
        dzm = cht.distance_int(Sign::neg(),Comp::k(),i,j,k,ts,ResistEval::yes,old);
        dtdzm = (phi[i][j][k]-ts)/dzm;
#ifndef VERSION_STABLE
        wmf = wpf;
#endif
      } else {
        if(dom->ibody().off(i,j,k-1)) {
          dtdzm = cht.gradt_ib(Sign::neg(),Comp::k(),i,j,k,old,phi);
        } else {
          dtdzm = (phi[i][j][k  ]-phi[i][j][k-1])/dzb(k);
        }
      }

      // dtdzp
      if(cht.interface(Sign::pos(),Comp::k(),i,j,k,old)) {
        real dzp, ts;
        dzp = cht.distance_int(Sign::pos(),Comp::k(),i,j,k,ts,ResistEval::yes,old);
        dtdzp = (ts-phi[i][j][k])/dzp;
#ifndef VERSION_STABLE
        wpf = wmf;
#endif
      } else {
        if(dom->ibody().off(i,j,k+1)) {
          dtdzp = cht.gradt_ib(Sign::pos(),Comp::k(),i,j,k,old,phi);
        } else {
          dtdzp = (phi[i][j][k+1]-phi[i][j][k  ])/dzt(k);
        }
      }
#endif

      uc  = 0.5*(umf+upf);
      vc  = 0.5*(vmf+vpf);
      wc  = 0.5*(wmf+wpf);

      udtdx = 0.5*(uc+fabs(uc))*dtdxm
            + 0.5*(uc-fabs(uc))*dtdxp;
      vdtdy = 0.5*(vc+fabs(vc))*dtdym
            + 0.5*(vc-fabs(vc))*dtdyp;
      wdtdz = 0.5*(wc+fabs(wc))*dtdzm
            + 0.5*(wc-fabs(wc))*dtdzp;

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
    if(cht.above_interface(i,j,k,old)) {
      c = cht.cpl(i,j,k);
    } else {
      c = cht.cpv(i,j,k);
    }
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }

}

