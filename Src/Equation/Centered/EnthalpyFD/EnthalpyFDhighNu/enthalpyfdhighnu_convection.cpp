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
  //setflag();

  /*-----------------------+
  |  finite volume method  |
  +-----------------------*/
  for_ijk(i,j,k) {
    //if(iflag[i][j][k]==0){
    //  continue;
    //}
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

    real phi_im2 = phi[i-2][j][k];
    real phi_im1 = phi[i-1][j][k];
    real phi_c   = phi[i  ][j][k];
    real phi_ip1 = phi[i+1][j][k];
    real phi_ip2 = phi[i+2][j][k];

    if( (clrold[i][j][k]-clrsurf)*(clrold[i-1][j][k]-clrsurf) < 0.0 ) {
      phi_im1 = phi_c;
      phi_im2 = phi_c;
    } else if( (clrold[i-1][j][k]-clrsurf)*(clrold[i-2][j][k]-clrsurf) < 0.0 )  {
      phi_im2 = phi_im1;
    }
    if( (clrold[i][j][k]-clrsurf)*(clrold[i+1][j][k]-clrsurf) < 0.0 ) {
      phi_ip1 = phi_c;
      phi_ip2 = phi_c;
    } else if( (clrold[i+1][j][k]-clrsurf)*(clrold[i+2][j][k]-clrsurf) < 0.0 ) {
      phi_ip2 = phi_ip1;
    }
    // flux at minus-surface, u is negative
    //phim = a_w*lim.limit(-umf, phi_ip,phi_c,phi_im);
    phim = a_w * ( lim.limit(+umf, phi_im2, phi_im1, phi_c)
                 + lim.limit(-umf, phi_ip1, phi_c, phi_im1));
    // flux at plus-surface, u is positive
    //phip = a_e*lim.limit(+upf, phi_im,phi_c,phi_ip);
    phip = a_e * ( lim.limit(+upf, phi_im1, phi_c, phi_ip1)
                 + lim.limit(-upf, phi_ip2, phi_ip1, phi_c));
 
    // comment out 2024.10.03 
    //if(i==si() && imin) phim = a_w*phi[i-1][j][k];
    //if(i==ei() && imax) phip = a_e*phi[i+1][j][k];

    (*conv)[i][j][k] += (umf*phim - upf*phip);
    //(*conv)[i+1][j][k] +=  upf*phip;
    //(*conv)[i-1][j][k] -=  umf*phim;

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

    real phi_jm2 = phi[i][j-2][k];
    real phi_jm1 = phi[i][j-1][k];
    real phi_c   = phi[i][j  ][k];
    real phi_jp1 = phi[i][j+1][k];
    real phi_jp2 = phi[i][j+2][k];

    if( (clrold[i][j][k]-clrsurf)*(clrold[i][j-1][k]-clrsurf) < 0.0 ) {
      phi_jm1 = phi_c;
      phi_jm2 = phi_c;
    } else if( (clrold[i][j-1][k]-clrsurf)*(clrold[i][j-2][k]-clrsurf) < 0.0 ) {
      phi_jm2 = phi_jm1;
    }
    if( (clrold[i][j][k]-clrsurf)*(clrold[i][j+1][k]-clrsurf) < 0.0 ) {
      phi_jp1 = phi_c;
      phi_jp2 = phi_c;
    } else if( (clrold[i][j][k]-clrsurf)*(clrold[i][j+1][k]-clrsurf) < 0.0 ) {
      phi_jp2 = phi_jp1;
    }

    //phim = a_s*lim.limit(-vmf, phi_jp,phi_c,phi_jp);
    phim = a_s * ( lim.limit(+vmf, phi_jm2, phi_jm1, phi_c)
                 + lim.limit(-vmf, phi_jp1, phi_c, phi_jm1));
    //phip = a_n*lim.limit(+vpf, phi_jm,phi_c,phi_jp);
    phip = a_n * ( lim.limit(+vpf, phi_jm1, phi_c, phi_jp1)
                 + lim.limit(-vpf, phi_jp2, phi_jp1, phi_c));

    //if(j==sj() && jmin) phim = a_s*phi[i][j-1][k];
    //if(j==ej() && jmax) phip = a_n*phi[i][j+1][k];

    (*conv)[i][j][k] += (vmf*phim - vpf*phip);
    //(*conv)[i][j+1][k] +=  vpf*phip;
    //(*conv)[i][j-1][k] -=  vmf*phim;

    } 
    { /////////
      //  w  //
      /////////
    wmf = (*u)[Comp::w()][i][j][k];    // w @ k-plus
    wpf = (*u)[Comp::w()][i][j][k+1];  // w @ k-minus

    real a_b = dSz(i,j,k);
    real a_t = dSz(i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_b *= dom->ibody().fSb(i,j,k);
      a_t *= dom->ibody().fSt(i,j,k);
    }

    real phi_km2 = phi[i][j][k-2];
    real phi_km1 = phi[i][j][k-1];
    real phi_c   = phi[i][j][k];
    real phi_kp1 = phi[i][j][k+1];
    real phi_kp2 = phi[i][j][k+2];

    if( (clrold[i][j][k]-clrsurf)*(clrold[i][j][k-1]-clrsurf) < 0.0 ) {
      phi_km1 = phi_c;
      phi_km2 = phi_c;
    } else if( (clrold[i][j][k-1]-clrsurf)*(clrold[i][j][k-2]-clrsurf) < 0.0 ) {
      phi_km2 = phi_km1;
    }
    if( (clrold[i][j][k]-clrsurf)*(clrold[i][j][k+1]-clrsurf) < 0.0 ) {
      phi_kp1 = phi_c;
      phi_kp2 = phi_c;
    } else if( (clrold[i][j][k+1]-clrsurf)*(clrold[i][j][k+2]-clrsurf) < 0.0 ) {
      phi_kp2 = phi_kp1;
    }

    //phim = a_b*lim.limit(-wmf, phi_kp, phi_c, phi_km);
    phim = a_b * ( lim.limit(+wmf, phi_km2, phi_km1, phi_c)
                 + lim.limit(-wmf, phi_kp1, phi_c, phi_km1));
    //phip = a_t*lim.limit(+wpf, phi_km, phi_c, phi_kp);
    phip = a_t * ( lim.limit(+wpf, phi_km1, phi_c, phi_kp1)
                 + lim.limit(-wpf, phi_kp2, phi_kp1, phi_c));

    //if(k==sk() && kmin) phim = a_b*phi[i][j][k-1];
    //if(k==ek() && kmax) phip = a_t*phi[i][j][k+1];

    (*conv)[i][j][k] += (wmf*phim - wpf*phip);
    //(*conv)[i][j][k+1] +=  wpf*phip;
    //(*conv)[i][j][k-1] -=  wmf*phim;
    
    }
  }

#if 0
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
#endif

  for_ijk(i,j,k) {
    //if(iflag[i][j][k]==0){
    //  continue;
    //}
    real divu = - dSx(i,j,k)*(*u)[Comp::u()][i]  [j]  [k]
                + dSx(i,j,k)*(*u)[Comp::u()][i+1][j]  [k]
                - dSy(i,j,k)*(*u)[Comp::v()][i]  [j]  [k]
                + dSy(i,j,k)*(*u)[Comp::v()][i]  [j+1][k]
                - dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k]
                + dSz(i,j,k)*(*u)[Comp::w()][i]  [j]  [k+1];
    (*conv)[i][j][k] += phi[i][j][k] * divu;
  }

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
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }
#endif

}
