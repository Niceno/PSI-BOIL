#include "enthalpyfdhighnu.h"
using namespace std;

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void EnthalpyFDhighNu::convection() {
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
void EnthalpyFDhighNu::convection(Scalar * conv) {

  phi.exchange();

  real phim, phip;
  real umf, upf, vmf, vpf, wmf, wpf;

  for_aijk(i,j,k)
    (*conv)[i][j][k] = 0.0;
#if 1
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

  //boil::plot->plot(*clr, iflag,"clr-phi-iflag", time->current_step());

  /* set flag */
  setflag();

  /*-----------------------+
  |  finite volume method  |
  +-----------------------*/
  for_ijk(i,j,k) {
    if(iflag[i][j][k]==0){
      continue;
    }
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

#if 0
    if(fabs(phi[i][j][k+1]-phi[i][j][k-1])>0.1){
      std::cout<<"k+1:k-1: "<<phi[i][j][k+1]<<" "<<phi[i][j][k-1]<<" "
               <<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      exit(0);
    }
#endif

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
    if(iflag[i][j][k]==0){
      continue;
    }
    real divu = - dSx(i,j,k)*(*u)[Comp::u()][i]  [j]  [k]
                + dSx(i,j,k)*(*u)[Comp::u()][i+1][j]  [k]
                - dSy(i,j,k)*(*u)[Comp::v()][i]  [j]  [k]
                + dSy(i,j,k)*(*u)[Comp::v()][i]  [j+1][k]
                - dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k]
                + dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k+1];
    (*conv)[i][j][k] += phi[i][j][k] * divu;
  }

#if 0
  if(boil::cart.iam()==0)
    std::cout<<"conv(32,34,34,0)= "<<(*conv)[32][34][34]<<"\n";
#endif

#if 1
  /*---------------------------+
  |  finite difference method  |
  +---------------------------*/
  for_ijk(i,j,k) {

    if(iflag[i][j][k]==0){

      real umf, upf, vmf, vpf, wmf, wpf, uc, vc, wc;
      real dtdxm, dtdxp, dtdym, dtdyp, dtdzm, dtdzp;
      real udtdx, vdtdy, wdtdz;

      //real clrc = (*clr)[i][j][k]; //shono edit
      real clrc = clrold[i][j][k];

      /* set velocity */
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
      //if((clrc-clrsurf)*((*clr)[i-1][j][k]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i-1][j][k]-clrsurf)<0.0){
        dtdxm = 0; //shono edit
      } else {
        dtdxm = (phi[i  ][j][k]-phi[i-1][j][k])/dxw(i);
      }

      // dtdxp
      //if((clrc-clrsurf)*((*clr)[i+1][j][k]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i+1][j][k]-clrsurf)<0.0){
        dtdxp = 0; //shono edit
      } else {
        dtdxp = (phi[i+1][j][k]-phi[i  ][j][k])/dxe(i);
      }

      // dtdym
      //if((clrc-clrsurf)*((*clr)[i][j-1][k]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i][j-1][k]-clrsurf)<0.0){
        dtdym = 0; //shono edit
      } else {
        dtdym = (phi[i][j  ][k]-phi[i][j-1][k])/dys(j);
      }

      // dtdyp
      //if((clrc-clrsurf)*((*clr)[i][j+1][k]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i][j+1][k]-clrsurf)<0.0){
        dtdyp = 0; //shono edit
      } else {
        dtdyp = (phi[i][j+1][k]-phi[i][j  ][k])/dyn(j);
      }

      // dtdzm
      //if((clrc-clrsurf)*((*clr)[i][j][k-1]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i][j][k-1]-clrsurf)<0.0){
        dtdzm = 0; //shono edit
      } else {
        dtdzm = (phi[i][j][k  ]-phi[i][j][k-1])/dzb(k);
      }

      // dtdzp
      //if((clrc-clrsurf)*((*clr)[i][j][k+1]-clrsurf)<0.0){ //shono edit
      if((clrc-clrsurf)*(clrold[i][j][k+1]-clrsurf)<0.0){ //shono edit
        dtdzp = 0; //shono edit
      } else {
        dtdzp = (phi[i][j][k+1]-phi[i][j][k  ])/dzt(k);
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

#if 0
  //if(boil::cart.iam()==0&&i==32&&j==34&&k==34)
  if(boil::cart.iam()==0&&i==32&&j==34&&k==34) {
    std::cout<<"FDM:conv(32,34,34,0)= "<<(*conv)[32][34][34]<<" "
             <<udtdx<<" "<<vdtdy<<" "<<wdtdz<<"\n";
    std::cout<<"if-state: "<<(clrc-clrsurf)*(clrold[i][j][k-1]-clrsurf)<<" "
             <<(clrc-clrsurf)*(clrold[i][j][k+1]-clrsurf)<<" "
             <<dtdzm<<" "<<dtdzp<<"\n";
    std::cout<<"phi: "<<phi[i][j][k-1]<<" "<<phi[i][j][k]<<" "<<phi[i][j][k+1]<<"\n";
    std::cout<<"clr: "<<clrold[i][j][k-1]<<" "<<clrold[i][j][k]<<" "<<clrold[i][j][k+1]<<"\n";
  }
#endif
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
    real c;
    if((*clr)[i][j][k]>=clrsurf){ //shono edit
    //if(clrold[i][j][k]>=clrsurf){
      c = cpl;
    } else {
      c = cpv;
    }
#if 0
    //std::cout<<"enthalpyfdhinu_convection\n";
    if (fabs((*conv)[i][j][k])>1e-16) {
      std::cout<<"convection: "<<(*conv)[i][j][k]<<" "<<i<<" "<<j<<" "<<k<<"\n";
      exit(0);
    }
#endif
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }
#endif


#if 0
  if(boil::cart.iam()==0)
    std::cout<<"End:conv(32,34,34,0)= "<<(*conv)[32][34][34]<<" "<<phi[32][34][34]<<" "<<iflag[32][34][34]<<"\n";
#endif
}
