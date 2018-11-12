#include "levelset.h"

/***************************************************************************//**
*  \brief Computes convection term using the last available velocities. 
*
*  \param conv - array into which convection term will be stored,
*
*  If called from new_time_step(), it will create \f$ \{C\}^{N-1} \f$; 
*  if called from inner iteration loop (from SIMPLE algorithm) it will compute
*  \f$ \{C\}^{N} \f$. 
*******************************************************************************/
void LevelSet::convection() {

  phi.exchange();

  real phim, phip;
  real umf, upf, vmf, vpf, wmf, wpf;

  for_aijk(i,j,k)
    stmp[i][j][k] = 0.0;

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

    stmp[i]  [j][k] += (umf*phim - upf*phip);
    stmp[i+1][j][k] +=  upf*phip;
    stmp[i-1][j][k] -=  umf*phim;
    }
    { /////////
      //     //
      //  v  //
      //     //
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

    stmp[i][j]  [k] += (vmf*phim - vpf*phip);
    stmp[i][j+1][k] +=  vpf*phip;
    stmp[i][j-1][k] -=  vmf*phim;
    } 
    { /////////
      //     //
      //  w  //
      //     //
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

    stmp[i][j][k]   += (wmf*phim - wpf*phip);
    stmp[i][j][k+1] +=  wpf*phip;
    stmp[i][j][k-1] -=  wmf*phim;
    }
  }

  for_ij(i,j) buff[i][j][ek()] = stmp[i][j][ek()+1]; 
  for_ij(i,j) buff[i][j][sk()] = stmp[i][j][sk()-1]; 
  buff.exchange(2);
  for_ij(i,j) stmp[i][j][ek()] += buff[i][j][ek()+1]; 
  for_ij(i,j) stmp[i][j][sk()] += buff[i][j][sk()-1]; 

  for_ik(i,k) buff[i][ej()][k] = stmp[i][ej()+1][k]; 
  for_ik(i,k) buff[i][sj()][k] = stmp[i][sj()-1][k]; 
  buff.exchange(1);
  for_ik(i,k) stmp[i][ej()][k] += buff[i][ej()+1][k]; 
  for_ik(i,k) stmp[i][sj()][k] += buff[i][sj()-1][k]; 

  for_jk(j,k) buff[ei()][j][k] = stmp[ei()+1][j][k]; 
  for_jk(j,k) buff[si()][j][k] = stmp[si()-1][j][k]; 
  buff.exchange(0);
  for_jk(j,k) stmp[ei()][j][k] += buff[ei()+1][j][k]; 
  for_jk(j,k) stmp[si()][j][k] += buff[si()-1][j][k];

#if 1
  for_ijk(i,j,k) {
    real a_w = dSx(i,j,k);
    real a_e = dSx(i,j,k);
    real a_s = dSy(i,j,k);
    real a_n = dSy(i,j,k);
    real a_b = dSz(i,j,k);
    real a_t = dSz(i,j,k);
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
    stmp[i][j][k] += phi[i][j][k] * divu;
  }
#endif

  /* update distance function */
  for_ijk(i,j,k){
    phi[i][j][k] += time->dt()*stmp[i][j][k]/dV(i,j,k);
  }

  insert_bc_dist(phi);
  phi.exchange_all();
}
