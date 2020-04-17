#include "topology.h"

/******************************************************************************/
void Topology::extrapolate(Scalar & sca, const Sign iext, 
                           const std::set<int> & testset) {
  extrapolate(sca,iext,testset,*iflag);
}

/******************************************************************************/
void Topology::extrapolate(Scalar & sca, const Sign iext,
                           const std::set<int> & testset, 
                           const ScalarInt & eflag) {
/***************************************************************************//**
*  \brief advect scalar variable sca in normal direction.
*         iext > 0 for extrapolate from vapor to liquid
*         iext < 0 for extrapolate from liquid to vapor
*         output : sca 
*         temporary: stmp, delta, stmp2
*******************************************************************************/
  /*----------------------------------------------------------+
  |  set flag                                                 |
  |  cell which will be extrapolated:                 stmp=0  |
  |  cell which will not be extrapolated (constant):  stmp=1  |
  +----------------------------------------------------------*/
  for_avijk(sca,i,j,k){
    int flagval = eflag[i][j][k];
    /* testset defines, which values of eflag mark cells to be extrapolated,
       the if-condition below is STL-standard way to check for being in set,
       equivalent to C++20 testset.contains(flagval) */
    if(testset.find(flagval) != testset.end()) {
      stmp[i][j][k] = 0; 
    /* otherwise */
    } else {
      stmp[i][j][k] = 1;
    }
    /* solid is always fixed */
    if(sca.domain()->ibody().off(i,j,k))
      stmp[i][j][k] = 1;
  }

#if 0
  if(iext<0) {
    stmp.exchange();
    boil::plot->plot((*clr),stmp,"clr-stmp",time->current_step());
    exit(0);
  }
#endif

  /*-------------------------------------------+
  |  advect sca                                |
  |  time scheme: euler implicit               |
  |  advection scheme: 1st order upwind        |
  |  iteration scheme: symmetric Gauss-Seidel  |
  +-------------------------------------------*/
  const real dtau = sca.domain()->dxyz_min();
  real isgn = real(iext);
  bool converged(false);

  // extrapolate from stmp=1 to stmp=0
  /* pseudo-time loop */
  for(int mstep=1; mstep<=mmax_ext; mstep++){
    /* right hand side */
    for_vijk(sca,i,j,k) {
      if(stmp[i][j][k]==0){

        real nxc = isgn*(*nx)[i][j][k];
        real nyc = isgn*(*ny)[i][j][k];
        real nzc = isgn*(*nz)[i][j][k];

        stmp2[i][j][k]  =-(0.5*(nxc+fabs(nxc))
                          *(sca[i][j][k]-sca[i-1][j][k])/sca.dxw(i));
        stmp2[i][j][k] -= (0.5*(nxc-fabs(nxc))
                          *(sca[i+1][j][k]-sca[i][j][k])/sca.dxe(i));
        stmp2[i][j][k] -= (0.5*(nyc+fabs(nyc))
                          *(sca[i][j][k]-sca[i][j-1][k])/sca.dys(j));
        stmp2[i][j][k] -= (0.5*(nyc-fabs(nyc))
                          *(sca[i][j+1][k]-sca[i][j][k])/sca.dyn(j));
        stmp2[i][j][k] -= (0.5*(nzc+fabs(nzc))
                          *(sca[i][j][k]-sca[i][j][k-1])/sca.dzb(k));
        stmp2[i][j][k] -= (0.5*(nzc-fabs(nzc))
                          *(sca[i][j][k+1]-sca[i][j][k])/sca.dzt(k));
      }
    }
    for_avijk(sca,i,j,k)
      delta[i][j][k]=0.0;

    /* loop for implicit scheme */
    for(int it=1; it<=4; it++){
      int ist,ied,iinc;
      int jst,jed,jinc;
      int kst,ked,kinc;
#if 1
      if(it%2==0){ist=sca.si();ied=sca.ei();iinc=1;}else{ist=sca.ei();ied=sca.si();iinc=-1;}
      if(it%2==0){jst=sca.sj();jed=sca.ej();jinc=1;}else{jst=sca.ej();jed=sca.sj();jinc=-1;}
      if(it%2==0){kst=sca.sk();ked=sca.ek();kinc=1;}else{kst=sca.ek();ked=sca.sk();kinc=-1;}
#else
      if(it%2==1){ist=sca.si();ied=sca.ei();iinc=1;}else{ist=sca.ei();ied=sca.si();iinc=-1;}
      if(it%2==1){jst=sca.sj();jed=sca.ej();jinc=1;}else{jst=sca.ej();jed=sca.sj();jinc=-1;}
      if(it%2==1){kst=sca.sk();ked=sca.ek();kinc=1;}else{kst=sca.ek();ked=sca.sk();kinc=-1;}
#endif
      for(int i=ist; i<=ied; i+=iinc){
      for(int j=jst; j<=jed; j+=jinc){
      for(int k=kst; k<=ked; k+=kinc){
        if(stmp[i][j][k]==0){
          real diag,drhs;
          real nxc = isgn*(*nx)[i][j][k];
          real nyc = isgn*(*ny)[i][j][k];
          real nzc = isgn*(*nz)[i][j][k];

          diag =1.0/dtau;
          diag+= 0.5*(nxc+fabs(nxc))/sca.dxw(i);
          drhs = 0.5*(nxc+fabs(nxc))*delta[i-1][j][k]/sca.dxw(i);

          diag-= 0.5*(nxc-fabs(nxc))/sca.dxe(i);
          drhs-= 0.5*(nxc-fabs(nxc))*delta[i+1][j][k]/sca.dxe(i);

          diag+= 0.5*(nyc+fabs(nyc))/sca.dys(j);
          drhs+= 0.5*(nyc+fabs(nyc))*delta[i][j-1][k]/sca.dys(j);

          diag-= 0.5*(nyc-fabs(nyc))/sca.dyn(j);
          drhs-= 0.5*(nyc-fabs(nyc))*delta[i][j+1][k]/sca.dyn(j);

          diag+= 0.5*(nzc+fabs(nzc))/sca.dzb(k);
          drhs+= 0.5*(nzc+fabs(nzc))*delta[i][j][k-1]/sca.dzb(k);

          diag-= 0.5*(nzc-fabs(nzc))/sca.dzt(k);
          drhs-= 0.5*(nzc-fabs(nzc))*delta[i][j][k+1]/sca.dzt(k);

          delta[i][j][k] = (stmp2[i][j][k]+drhs)/diag;
        }
      }}}
    }

    /* update sca */
    real errnorm(0.0);
    for_vijk(sca,i,j,k) {
      if(stmp[i][j][k]==0){
        real diff = delta[i][j][k];
        real & modd = sca[i][j][k];
        real err = fabs(diff/(modd+boil::pico));
        if(err>errnorm)
          errnorm = err;
        modd += diff;
      }
    }
    sca.exchange();
    boil::cart.max_real(&errnorm);
    //boil::oout<<mstep<<" "<<errnorm<<boil::endl;
  
    if(errnorm<tol_ext) {
      converged = true;
      //boil::oout<<"Topology::ext_grad converged after "<<mstep<<" steps, final rel. error: "<<errnorm<<"\n";
      break;
    }

  }

  if(!converged)
    boil::oout<<"Topology::ext_grad did not converge after "<<mmax_ext<<" steps! \n";

  return;
}

