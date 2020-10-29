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

  /*-----------------------+
  |  finite volume method  |
  +-----------------------*/
  for_ijk(i,j,k) {
    { /////////
      //  u  //
      /////////
#ifdef USE_PHASIC_VELOCITIES
  #ifdef CNEW
    if(cht.topo->above_interface(i,j,k)) {
  #else
    if(cht.topo->above_interface_old(i,j,k)) {
  #endif
      umf = (*uliq)[Comp::u()][i]  [j][k];  // u @ imin
      upf = (*uliq)[Comp::u()][i+1][j][k];  // u @ imax
    } else {
      umf = (*ugas)[Comp::u()][i]  [j][k];  // u @ imin
      upf = (*ugas)[Comp::u()][i+1][j][k];  // u @ imax
    }
#else
    umf = (*u)[Comp::u()][i]  [j][k];  // u @ imin
    upf = (*u)[Comp::u()][i+1][j][k];  // u @ imax
#endif
    
    real a_w = dSx(Sign::neg(),i,j,k);
    real a_e = dSx(Sign::pos(),i,j,k);
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
#ifdef USE_PHASIC_VELOCITIES
  #ifdef CNEW
    if(cht.topo->above_interface(i,j,k)) {
  #else
    if(cht.topo->above_interface_old(i,j,k)) {
  #endif
      vmf = (*uliq)[Comp::v()][i][j]  [k];  // v @ jmin
      vpf = (*uliq)[Comp::v()][i][j+1][k];  // v @ jmax
    } else {
      vmf = (*ugas)[Comp::v()][i][j]  [k];  // v @ jmin
      vpf = (*ugas)[Comp::v()][i][j+1][k];  // v @ jmax
    }
#else
    vmf = (*u)[Comp::v()][i][j]  [k];  // v @ jmin
    vpf = (*u)[Comp::v()][i][j+1][k];  // v @ jmax
#endif

    real a_s = dSy(Sign::neg(),i,j,k);
    real a_n = dSy(Sign::pos(),i,j,k);
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
#ifdef USE_PHASIC_VELOCITIES
  #ifdef CNEW
    if(cht.topo->above_interface(i,j,k)) {
  #else
    if(cht.topo->above_interface_old(i,j,k)) {
  #endif
      wmf = (*uliq)[Comp::w()][i][j][k];    // w @ kmin
      wpf = (*uliq)[Comp::w()][i][j][k+1];  // w @ kmax
    } else {
      wmf = (*ugas)[Comp::w()][i][j][k];    // w @ kmin
      wpf = (*ugas)[Comp::w()][i][j][k+1];  // w @ kmax
    }
#else
    wmf = (*u)[Comp::w()][i][j][k];    // w @ kmin
    wpf = (*u)[Comp::w()][i][j][k+1];  // w @ kmax
#endif

    real a_b = dSz(Sign::neg(),i,j,k);
    real a_t = dSz(Sign::pos(),i,j,k);
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
#if 0
    real divu = - dSx(Sign::neg(),i,j,k)*(*u)[Comp::u()][i]  [j]  [k]
                + dSx(Sign::pos(),i,j,k)*(*u)[Comp::u()][i+1][j]  [k]
                - dSy(Sign::neg(),i,j,k)*(*u)[Comp::v()][i]  [j]  [k]
                + dSy(Sign::pos(),i,j,k)*(*u)[Comp::v()][i]  [j+1][k]
                - dSz(Sign::neg(),i,j,k)*(*u)[Comp::w()][i]  [j]  [k]
                + dSz(Sign::pos(),i,j,k)*(*u)[Comp::w()][i]  [j]  [k+1];
#else
  #ifdef USE_PHASIC_VELOCITIES
    real divu;
  #ifdef CNEW
    if(cht.topo->above_interface(i,j,k)) {
  #else
    if(cht.topo->above_interface_old(i,j,k)) {
  #endif
      divu = uliq->outflow(i,j,k);
    } else {
      divu = ugas->outflow(i,j,k);
    }
  #else
    real divu = u->outflow(i,j,k);
  #endif
#endif
    (*conv)[i][j][k] += phi[i][j][k] * divu;
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

#ifdef USE_PHASIC_VELOCITIES
  #ifdef CNEW
      if(cht.topo->above_interface(i,j,k)) {
  #else
      if(cht.topo->above_interface_old(i,j,k)) {
  #endif
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
#else
      // u
      umf = (*u)[Comp::u()][i]  [j][k];
      upf = (*u)[Comp::u()][i+1][j][k];

      // v
      vmf = (*u)[Comp::v()][i][j]  [k];
      vpf = (*u)[Comp::v()][i][j+1][k];

      // w
      wmf = (*u)[Comp::w()][i][j][k];
      wpf = (*u)[Comp::w()][i][j][k+1];
#endif

#if 1
  #ifdef CNEW
      dtdxm = dtdxp = cht.gradt1D(false,Comp::i(),i,j,k,ao_conv,false,Old::no);
      dtdym = dtdyp = cht.gradt1D(false,Comp::j(),i,j,k,ao_conv,false,Old::no);
      dtdzm = dtdzp = cht.gradt1D(false,Comp::k(),i,j,k,ao_conv,false,Old::no);
  #else
      dtdxm = dtdxp = cht.gradt1D(false,Comp::i(),i,j,k,ao_conv,false,Old::yes);
      dtdym = dtdyp = cht.gradt1D(false,Comp::j(),i,j,k,ao_conv,false,Old::yes);
      dtdzm = dtdzp = cht.gradt1D(false,Comp::k(),i,j,k,ao_conv,false,Old::yes);
  #endif
#else
      // dtdxm
  #ifdef CNEW
      if(cht.interface(Sign::neg(),Comp::i(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::neg(),Comp::i(),i,j,k)) {
  #endif
        real dxm, ts;
  #ifdef CNEW
        dxm = cht.distance_int_x(Sign::neg(),i,j,k,ts);
  #else
        dxm = cht.distance_int_x_old(Sign::neg(),i,j,k,ts);
  #endif      
        dtdxm = (phi[i][j][k]-ts)/dxm;
#ifndef VERSION_STABLE
        umf = upf;
#endif
      } else {
        if(dom->ibody().off(i-1,j,k)) {
          dtdxm = cht.gradt_ib(Sign::neg(),Comp::i(),i,j,k,Old::yes,phi);
        } else {
          dtdxm = (phi[i  ][j][k]-phi[i-1][j][k])/dxw(i);
        }
      }

      // dtdxp
  #ifdef CNEW
      if(cht.interface(Sign::pos(),Comp::i(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::pos(),Comp::i(),i,j,k)) {
  #endif
        real dxp, ts;
  #ifdef CNEW
        dxp = cht.distance_int_x(Sign::pos(),i,j,k,ts);
  #else
        dxp = cht.distance_int_x_old(Sign::pos(),i,j,k,ts);
  #endif      
        dtdxp = (ts-phi[i][j][k])/dxp;
#ifndef VERSION_STABLE
        upf = umf;
#endif
      } else {
        if(dom->ibody().off(i+1,j,k)) {
          dtdxp = cht.gradt_ib(Sign::pos(),Comp::i(),i,j,k,Old::yes,phi);
        } else {
          dtdxp = (phi[i+1][j][k]-phi[i  ][j][k])/dxe(i);
        }
      }

      // dtdym
  #ifdef CNEW
      if(cht.interface(Sign::neg(),Comp::j(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::neg(),Comp::j(),i,j,k)) {
  #endif
        real dym, ts;
  #ifdef CNEW
        dym = cht.distance_int_y(Sign::neg(),i,j,k,ts);
  #else
        dym = cht.distance_int_y_old(Sign::neg(),i,j,k,ts);
  #endif      
        dtdym = (phi[i][j][k]-ts)/dym;
#ifndef VERSION_STABLE
        vmf = vpf;
#endif
      } else {
        if(dom->ibody().off(i,j-1,k)) {
          dtdym = cht.gradt_ib(Sign::neg(),Comp::j(),i,j,k,Old::yes,phi);
        } else {
          dtdym = (phi[i][j  ][k]-phi[i][j-1][k])/dys(j);
        }
      }

      // dtdyp
  #ifdef CNEW
      if(cht.interface(Sign::pos(),Comp::j(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::pos(),Comp::j(),i,j,k)) {
  #endif
        real dyp, ts;
  #ifdef CNEW
        dyp = cht.distance_int_y(Sign::pos(),i,j,k,ts);
  #else
        dyp = cht.distance_int_y_old(Sign::pos(),i,j,k,ts);
  #endif      
        dtdyp = (ts-phi[i][j][k])/dyp;
#ifndef VERSION_STABLE
        vpf = vmf;
#endif
      } else {
        if(dom->ibody().off(i,j+1,k)) {
          dtdyp = cht.gradt_ib(Sign::pos(),Comp::j(),i,j,k,Old::yes,phi);
        } else {
          dtdyp = (phi[i][j+1][k]-phi[i][j  ][k])/dyn(j);
        }
      }

      // dtdzm
  #ifdef CNEW
      if(cht.interface(Sign::neg(),Comp::k(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::neg(),Comp::k(),i,j,k)) {
  #endif
        real dzm, ts;
  #ifdef CNEW
        dzm = cht.distance_int_z(Sign::neg(),i,j,k,ts);
  #else
        dzm = cht.distance_int_z_old(Sign::neg(),i,j,k,ts);
  #endif      
        dtdzm = (phi[i][j][k]-ts)/dzm;
#ifndef VERSION_STABLE
        wmf = wpf;
#endif
      } else {
        if(dom->ibody().off(i,j,k-1)) {
          dtdzm = cht.gradt_ib(Sign::neg(),Comp::k(),i,j,k,Old::yes,phi);
        } else {
          dtdzm = (phi[i][j][k  ]-phi[i][j][k-1])/dzb(k);
        }
      }

      // dtdzp
  #ifdef CNEW
      if(cht.interface(Sign::pos(),Comp::k(),i,j,k)) {
  #else
      if(cht.interface_old(Sign::pos(),Comp::k(),i,j,k)) {
  #endif
        real dzp, ts;
  #ifdef CNEW
        dzp = cht.distance_int_z(Sign::pos(),i,j,k,ts);
  #else
        dzp = cht.distance_int_z_old(Sign::pos(),i,j,k,ts);
  #endif      
        dtdzp = (ts-phi[i][j][k])/dzp;
#ifndef VERSION_STABLE
        wpf = wmf;
#endif
      } else {
        if(dom->ibody().off(i,j,k+1)) {
          dtdzp = cht.gradt_ib(Sign::pos(),Comp::k(),i,j,k,Old::yes,phi);
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
  #ifdef CNEW
    if(cht.topo->above_interface(i,j,k)) {
  #else
    if(cht.topo->above_interface_old(i,j,k)) {
  #endif
      c = cht.cpl(i,j,k);
    } else {
      c = cht.cpv(i,j,k);
    }
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }

}

