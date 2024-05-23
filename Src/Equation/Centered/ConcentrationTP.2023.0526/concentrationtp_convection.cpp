#include "concentrationtp.h"

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void ConcentrationTP::convection() {
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
void ConcentrationTP::convection(Scalar * conv) {

  phi.exchange();

  real phim, phip, flxm, flxp, gfxm, gfxp;
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

  /*-----------------------+
  |  finite volume method  |
  +-----------------------*/
  for_ijk(i,j,k) {
    { /////////
      //     //
      //  u  //
      //     //
      /////////
    Comp mcomp = Comp::u();

    umf = (*u)[mcomp][i  ][j][k];  /* u @ imin */
    upf = (*u)[mcomp][i+1][j][k];  /* u @ imax */

    if(matter_sig==Sign::neg()) {
      //std::cout<<"convection: "<<i<<" "<<j<<" "<<k<<"\n";
      //exit(0);
      /* liquid flux */
      flxm = (*colorflow)[mcomp][i  ][j][k]*time->dti();
      flxp = (*colorflow)[mcomp][i+1][j][k]*time->dti();

      /* total volume flux */
      gfxm = dSx(Sign::neg(),i,j,k)*umf;
      gfxp = dSx(Sign::pos(),i,j,k)*upf;

      if(dom->ibody().cut(i,j,k)) {
        gfxm *= dom->ibody().fSw(i,j,k);
        gfxp *= dom->ibody().fSe(i,j,k);
      }
      /* gas flux */
      gfxm = gfxm - flxm;
      gfxp = gfxp - flxp;
    } else {
      //std::cout<<"convection: "<<i<<" "<<j<<" "<<k<<"\n";
      //exit(0);
      gfxm = (*colorflow)[mcomp][i  ][j][k]*time->dti();
      gfxp = (*colorflow)[mcomp][i+1][j][k]*time->dti();
    }
 
    phim = lim.limit(-gfxm, (phi)[i+1][j][k],(phi)[i][j][k],(phi)[i-1][j][k]);
    if(i==si() && imin) phim = (phi)[i-1][j][k];
    phip = lim.limit(+gfxp, (phi)[i-1][j][k],(phi)[i][j][k],(phi)[i+1][j][k]);
    if(i==ei() && imax) phip = (phi)[i+1][j][k];
#if 0
    if(phim>1.01||phip>1.01){
      std::cout<<"convection:phim,phip "<<phim<<" "<<phip<<"\n";
      std::cout<<(phi)[i+1][j][k]<<" "<<(phi)[i][j][k]<<" "<<(phi)[i-1][j][k]<<"\n";
      exit(0);
    }
#endif

    /* transport in gas */
    (*conv)[i  ][j][k] += (gfxm*phim - gfxp*phip);
    (*conv)[i+1][j][k] += gfxp*phip;
    (*conv)[i-1][j][k] -= gfxm*phim;

    }
    { /////////
      //     //
      //  v  //
      //     //
      /////////
    Comp mcomp = Comp::v();

    vmf = (*u)[mcomp][i][j  ][k];  /* v @ jmin */
    vpf = (*u)[mcomp][i][j+1][k];  /* v @ jmax */

    if(matter_sig==Sign::neg()) {
      /* liquid flux */
      flxm = (*colorflow)[mcomp][i][j  ][k]*time->dti();
      flxp = (*colorflow)[mcomp][i][j+1][k]*time->dti();

      /* total volume flux */
      gfxm = dSy(Sign::neg(),i,j,k)*vmf;
      gfxp = dSy(Sign::pos(),i,j,k)*vpf;
 
      if(dom->ibody().cut(i,j,k)) {
        gfxm *= dom->ibody().fSs(i,j,k);
        gfxp *= dom->ibody().fSn(i,j,k);
      }
  
      /* gas flux */
      gfxm = gfxm - flxm;
      gfxp = gfxp - flxp;
    } else {
      gfxm = (*colorflow)[mcomp][i][j  ][k]*time->dti();
      gfxp = (*colorflow)[mcomp][i][j+1][k]*time->dti();
    }

    phim = lim.limit(-gfxm, (phi)[i][j+1][k],(phi)[i][j][k],(phi)[i][j-1][k]);
    if(j==sj() && jmin) phim = (phi)[i][j-1][k];   
    phip = lim.limit(+gfxp, (phi)[i][j-1][k],(phi)[i][j][k],(phi)[i][j+1][k]);
    if(j==ej() && jmax) phip = (phi)[i][j+1][k];

    /* transport in gas */
    (*conv)[i][j  ][k] += (gfxm*phim - gfxp*phip);
    (*conv)[i][j+1][k] += gfxp*phip;
    (*conv)[i][j-1][k] -= gfxm*phim;

    }

    { /////////
      //     //
      //  w  //
      //     //
      /////////
    Comp mcomp = Comp::w();

    wmf = (*u)[mcomp][i][j][k  ];  /* w @ kmin */
    wpf = (*u)[mcomp][i][j][k+1];  /* w @ kmax */

    if(matter_sig==Sign::neg()) {
      /* liquid flux */
      flxm = (*colorflow)[mcomp][i][j][k  ]*time->dti();
      flxp = (*colorflow)[mcomp][i][j][k+1]*time->dti();
  
      /* total volume flux */
      gfxm = dSz(Sign::neg(),i,j,k)*wmf;
      gfxp = dSz(Sign::pos(),i,j,k)*wpf;
 
      if(dom->ibody().cut(i,j,k)) {
        gfxm *= dom->ibody().fSb(i,j,k);
        gfxp *= dom->ibody().fSt(i,j,k);
      }

      /* gas flux */
      gfxm = gfxm - flxm;
      gfxp = gfxp - flxp;
    } else {
      gfxm = (*colorflow)[mcomp][i][j][k  ]*time->dti();
      gfxp = (*colorflow)[mcomp][i][j][k+1]*time->dti();
    }

    phim = lim.limit(-gfxm, (phi)[i][j][k+1], (phi)[i][j][k], (phi)[i][j][k-1]);
    if(k==sk() && kmin) phim = (phi)[i][j][k-1];
    phip = lim.limit(+gfxp, (phi)[i][j][k-1], (phi)[i][j][k], (phi)[i][j][k+1]);
    if(k==ek() && kmax) phip = (phi)[i][j][k+1];

    /* transport in gas */
    (*conv)[i][j][k  ] += (gfxm*phim - gfxp*phip);
    (*conv)[i][j][k+1] +=  gfxp*phip;
    (*conv)[i][j][k-1] -=  gfxm*phim;

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
    real r = rho_dif->value(i,j,k);
    (*conv)[i][j][k] = r * (*conv)[i][j][k];
  }

}
