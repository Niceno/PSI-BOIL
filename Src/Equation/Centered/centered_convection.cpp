#include "centered.h"

#define NEW_CONSERVATIVE_FORM true

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void Centered::convection() {
  convection(&cnew, flu->rho());
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
void Centered::convection(Scalar * conv, const Property * prop) {

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
        phi.bc().type(b)==BndType::insert() ) {
      if(phi.bc().direction(b) == Dir::imin() && 
         dom->coord(Comp::i()) == 0)   
        imin=true;
      if(phi.bc().direction(b) == Dir::imax() && 
         dom->coord(Comp::i()) == dom->dim(Comp::i())-1) 
        imax=true;
      if(phi.bc().direction(b) == Dir::jmin() && 
         dom->coord(Comp::j()) == 0)              
        jmin=true;
      if(phi.bc().direction(b) == Dir::jmax() && 
         dom->coord(Comp::j()) == dom->dim(Comp::j())-1) 
        jmax=true;
      if(phi.bc().direction(b) == Dir::kmin() && 
         dom->coord(Comp::k()) == 0)              
        kmin=true;
      if(phi.bc().direction(b) == Dir::kmax() && 
         dom->coord(Comp::k()) == dom->dim(Comp::k())-1) 
        kmax=true;
    }
  }

  /*--------------------+
  |  inside the domain  |
  +--------------------*/
  for_ijk(i,j,k) {
    { /////////
      //     //
      //  u  //
      //     //
      /////////
    umf = (*u)[Comp::u()][i]  [j][k];  // u @ imin
    upf = (*u)[Comp::u()][i+1][j][k];  // u @ imax
    
    real a_w = dSx(Sign::neg(),i,j,k);
    real a_e = dSx(Sign::pos(),i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_w *= dom->ibody().fSw(i,j,k);
      a_e *= dom->ibody().fSe(i,j,k);
    }

    #if !NEW_CONSERVATIVE_FORM
    const real rm = prop->value(Comp::u(),i  ,j,k);
    const real rp = prop->value(Comp::u(),i+1,j,k);
 
    phim = rm*a_w*lim.limit(-umf,phi[i+1][j][k],phi[i][j][k],phi[i-1][j][k]);
    phip = rp*a_e*lim.limit(+upf,phi[i-1][j][k],phi[i][j][k],phi[i+1][j][k]);

    if(i==si() && imin) phim = rm*a_w*phi[i-1][j][k];
    if(i==ei() && imax) phip = rp*a_e*phi[i+1][j][k];
    #else
    phim = a_w*lim.limit(-umf, phi[i+1][j][k],phi[i][j][k],phi[i-1][j][k]);
    phip = a_e*lim.limit(+upf, phi[i-1][j][k],phi[i][j][k],phi[i+1][j][k]);
  
    if(i==si() && imin) phim = a_w*phi[i-1][j][k];
    if(i==ei() && imax) phip = a_e*phi[i+1][j][k];
    #endif

    (*conv)[i]  [j][k] += (umf*phim - upf*phip);
    (*conv)[i+1][j][k] +=  upf*phip;
    (*conv)[i-1][j][k] -=  umf*phim;
    }
    { /////////
      //     //
      //  v  //
      //     //
      /////////
    vmf = (*u)[Comp::v()][i][j]  [k];  // v @ jmin
    vpf = (*u)[Comp::v()][i][j+1][k];  // v @ jmax

    real a_s = dSy(Sign::neg(),i,j,k);
    real a_n = dSy(Sign::pos(),i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_s *= dom->ibody().fSs(i,j,k);
      a_n *= dom->ibody().fSn(i,j,k);
    }

    #if !NEW_CONSERVATIVE_FORM
    const real rm = prop->value(Comp::u(),i,j  ,k);
    const real rp = prop->value(Comp::u(),i,j+1,k);
    phim = rm*a_s*lim.limit(-vmf,phi[i][j+1][k],phi[i][j][k],phi[i][j-1][k]);
    phip = rp*a_n*lim.limit(+vpf,phi[i][j-1][k],phi[i][j][k],phi[i][j+1][k]);

    if(j==sj() && jmin) phim = rm*a_s*phi[i][j-1][k];
    if(j==ej() && jmax) phip = rp*a_n*phi[i][j+1][k];
    #else
    phim = a_s*lim.limit(-vmf, phi[i][j+1][k],phi[i][j][k],phi[i][j-1][k]);
    phip = a_n*lim.limit(+vpf, phi[i][j-1][k],phi[i][j][k],phi[i][j+1][k]);

    if(j==sj() && jmin) phim = a_s*phi[i][j-1][k];
    if(j==ej() && jmax) phip = a_n*phi[i][j+1][k];
    #endif

    (*conv)[i][j]  [k] += (vmf*phim - vpf*phip);
    (*conv)[i][j+1][k] +=  vpf*phip;
    (*conv)[i][j-1][k] -=  vmf*phim;
    } 
    { /////////
      //     //
      //  w  //
      //     //
      /////////
    wmf = (*u)[Comp::w()][i][j][k];    // w @ kmin
    wpf = (*u)[Comp::w()][i][j][k+1];  // w @ kmax

    real a_b = dSz(Sign::neg(),i,j,k);
    real a_t = dSz(Sign::pos(),i,j,k);
    if(dom->ibody().cut(i,j,k)) {
      a_b *= dom->ibody().fSb(i,j,k);
      a_t *= dom->ibody().fSt(i,j,k);
    }

    #if !NEW_CONSERVATIVE_FORM
    const real rm = prop->value(Comp::u(),i,j,k  );
    const real rp = prop->value(Comp::u(),i,j,k+1);
    phim = rm*a_b*lim.limit(-wmf,phi[i][j][k+1],phi[i][j][k],phi[i][j][k-1]);
    phip = rp*a_t*lim.limit(+wpf,phi[i][j][k-1],phi[i][j][k],phi[i][j][k+1]);

    if(k==sk() && kmin) phim = rm*a_b*phi[i][j][k-1];
    if(k==ek() && kmax) phip = rp*a_t*phi[i][j][k+1];
    #else
    phim = a_b*lim.limit(-wmf, phi[i][j][k+1], phi[i][j][k], phi[i][j][k-1]);
    phip = a_t*lim.limit(+wpf, phi[i][j][k-1], phi[i][j][k], phi[i][j][k+1]);

    if(k==sk() && kmin) phim = a_b*phi[i][j][k-1];
    if(k==ek() && kmax) phip = a_t*phi[i][j][k+1];
    #endif

    (*conv)[i][j][k]   += (wmf*phim - wpf*phip);
    (*conv)[i][j][k+1] +=  wpf*phip;
    (*conv)[i][j][k-1] -=  wmf*phim;
    }
  }

  for_ij(i,j) buff[i][j][ek()] = (*conv)[i][j][ek()+1]; 
  for_ij(i,j) buff[i][j][sk()] = (*conv)[i][j][sk()-1]; 
  buff.exchange(2);
  for_ij(i,j) (*conv)[i][j][ek()] += buff[i][j][ek()+1]; 
  for_ij(i,j) (*conv)[i][j][sk()] += buff[i][j][sk()-1]; 

  for_ik(i,k) buff[i][ej()][k] = (*conv)[i][ej()+1][k]; 
  for_ik(i,k) buff[i][sj()][k] = (*conv)[i][sj()-1][k]; 
  buff.exchange(1);
  for_ik(i,k) (*conv)[i][ej()][k] += buff[i][ej()+1][k]; 
  for_ik(i,k) (*conv)[i][sj()][k] += buff[i][sj()-1][k]; 

  for_jk(j,k) buff[ei()][j][k] = (*conv)[ei()+1][j][k]; 
  for_jk(j,k) buff[si()][j][k] = (*conv)[si()-1][j][k]; 
  buff.exchange(0);
  for_jk(j,k) (*conv)[ei()][j][k] += buff[ei()+1][j][k]; 
  for_jk(j,k) (*conv)[si()][j][k] += buff[si()-1][j][k]; 

  #if NEW_CONSERVATIVE_FORM
  for_ijk(i,j,k) {

    real a_w = dSx(Sign::neg(),i,j,k);
    real a_e = dSx(Sign::pos(),i,j,k);
    real a_s = dSy(Sign::neg(),i,j,k);
    real a_n = dSy(Sign::pos(),i,j,k);
    real a_b = dSz(Sign::neg(),i,j,k);
    real a_t = dSz(Sign::pos(),i,j,k);

    if(dom->ibody().cut(i,j,k)) {
      a_w *= dom->ibody().fSw(i,j,k);
      a_e *= dom->ibody().fSe(i,j,k);
      a_b *= dom->ibody().fSb(i,j,k);
      a_t *= dom->ibody().fSt(i,j,k);
      a_s *= dom->ibody().fSs(i,j,k);
      a_n *= dom->ibody().fSn(i,j,k);
    }

    real divu = - a_w * (*u)[Comp::u()][i]  [j]  [k]
                + a_e * (*u)[Comp::u()][i+1][j]  [k]
                - a_s * (*u)[Comp::v()][i]  [j]  [k]
                + a_n * (*u)[Comp::v()][i]  [j+1][k]
                - a_b * (*u)[Comp::w()][i]  [j]  [k]
                + a_t * (*u)[Comp::w()][i]  [j]  [k+1];
    (*conv)[i][j][k] += phi[i][j][k] * divu;
  }

  for_ijk(i,j,k) {
    (*conv)[i][j][k] = prop->value(i,j,k) * (*conv)[i][j][k];
  }
  #endif

  //if(time->current_step() % 100 == 0)
  //  boil::plot->plot(*conv, "conv", time->current_step());
}
