#include "phasechange.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange::ext_gradt(Scalar & sca, const int iext) {
/***************************************************************************//**
*  \brief advect scalar variable sca in normal direction.
*         iext =  1 for extrapolate from vapor to liquid
*         iext = -1 for extrapolate from liquid to vapor
*         output : sca 
*         tempolary: stmp2
*******************************************************************************/
  /*----------------------------------------------------------+
  |  set flag                                                 |
  |  cell which will be extrapolated:                 stmp=0  |
  |  cell which will not be extrapolated (constant):  stmp=1  |
  +----------------------------------------------------------*/
  if(iext==1){                 // extrapolate from vapor to liquid
    for_aijk(i,j,k){
      if( iflag[i][j][k] == 1 ){   // liquid cell next to interface
        stmp[i][j][k]=0;
      } else {                 // otherwise
        stmp[i][j][k]=1;
      }
      if(dom->ibody().off(i,j,k))stmp[i][j][k]=1; // fix value in solid
    }
  } else {                     // extrapolate from liqid to vapor
    for_aijk(i,j,k){
      if( iflag[i][j][k] == -1 ){  // vapor cell next to interface
        stmp[i][j][k]=0;
      } else {                 // otherwise
        stmp[i][j][k]=1;
      }
      if(dom->ibody().off(i,j,k))stmp[i][j][k]=1; // fix value in solid
    }
  }
  //stmp.exchange();
  //boil::plot->plot(clr,stmp,"clr-stmp",time->current_step());

  /*-------------------------------------------+
  |  advect sca                                |
  |  time scheme: euler implicit               |
  |  advection scheme: 1st order upwind        |
  |  iteration scheme: symmetric Gauss-Seidel  |
  +-------------------------------------------*/
  const real dtau = dxmin;
  real isgn = real(iext);
  int mmax=8;
  //int mmax=20;

  // extrapolate from stmp=0 to stmp=1
  /* pseudo-time loop */
  for(int mstep=1; mstep<=mmax; mstep++){
    /* right hand side */
    for_ijk(i,j,k) {
      if(int(stmp[i][j][k])==0){

        real nxc = isgn*nx[i][j][k];
        real nyc = isgn*ny[i][j][k];
        real nzc = isgn*nz[i][j][k];

        stmp2[i][j][k]  =-(0.5*(nxc+fabs(nxc))
                          *(sca[i][j][k]-sca[i-1][j][k])/dxw(i));
        stmp2[i][j][k] -= (0.5*(nxc-fabs(nxc))
                          *(sca[i+1][j][k]-sca[i][j][k])/dxe(i));
        stmp2[i][j][k] -= (0.5*(nyc+fabs(nyc))
                          *(sca[i][j][k]-sca[i][j-1][k])/dys(j));
        stmp2[i][j][k] -= (0.5*(nyc-fabs(nyc))
                          *(sca[i][j+1][k]-sca[i][j][k])/dyn(j));
        stmp2[i][j][k] -= (0.5*(nzc+fabs(nzc))
                          *(sca[i][j][k]-sca[i][j][k-1])/dzb(k));
        stmp2[i][j][k] -= (0.5*(nzc-fabs(nzc))
                          *(sca[i][j][k+1]-sca[i][j][k])/dzt(k));
      }
    }
    for_aijk(i,j,k)
      delta[i][j][k]=0.0;

    /* loop for implicit scheme */
    for(int it=1; it<=4; it++){
      int ist,ied,iinc;
      int jst,jed,jinc;
      int kst,ked,kinc;
      if(it%2==0){ist=si();ied=ei();iinc=1;}else{ist=ei();ied=si();iinc=-1;}
      if(it%2==0){jst=sj();jed=ej();jinc=1;}else{jst=ej();jed=sj();jinc=-1;}
      if(it%2==0){kst=sk();ked=ek();kinc=1;}else{kst=ek();ked=sk();kinc=-1;}
      for(int i=ist; i<=ied; i+=iinc){
      for(int j=jst; j<=jed; j+=jinc){
      for(int k=kst; k<=ked; k+=kinc){
        if(int(stmp[i][j][k])==0){
          real diag,drhs;
          real nxc = isgn*nx[i][j][k];
          real nyc = isgn*ny[i][j][k];
          real nzc = isgn*nz[i][j][k];

          diag =1.0/dtau;
          diag+= 0.5*(nxc+fabs(nxc))/dxw(i);
          drhs = 0.5*(nxc+fabs(nxc))*delta[i-1][j][k]/dxw(i);

          diag-= 0.5*(nxc-fabs(nxc))/dxe(i);
          drhs-= 0.5*(nxc-fabs(nxc))*delta[i+1][j][k]/dxe(i);

          diag+= 0.5*(nyc+fabs(nyc))/dys(j);
          drhs+= 0.5*(nyc+fabs(nyc))*delta[i][j-1][k]/dys(j);

          diag-= 0.5*(nyc-fabs(nyc))/dyn(j);
          drhs-= 0.5*(nyc-fabs(nyc))*delta[i][j+1][k]/dyn(j);

          diag+= 0.5*(nzc+fabs(nzc))/dzb(k);
          drhs+= 0.5*(nzc+fabs(nzc))*delta[i][j][k-1]/dzb(k);

          diag-= 0.5*(nzc-fabs(nzc))/dzt(k);
          drhs-= 0.5*(nzc-fabs(nzc))*delta[i][j][k+1]/dzt(k);

          delta[i][j][k] = (stmp2[i][j][k]+drhs)/diag;
        }
      }}}
    }

    /* update sca */
    for_ijk(i,j,k) {
      if(int(stmp[i][j][k])==0){
        sca[i][j][k] += delta[i][j][k];
      }
    }
    sca.exchange();
  }

  return;
}
