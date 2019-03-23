#include "vof.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void VOF::ext_vel(Scalar & sca, const Scalar & eflag, const int iext) {
/***************************************************************************//**
*  \brief advect scalar variable sca in normal direction.
*         iext =  1 for extrapolate from vapor to liquid
*         iext = -1 for extrapolate from liquid to vapor
*         output : sca 
*         temporary: stmp2, stmp3
*******************************************************************************/
  /*-------------------------------------------------------------+
  |  cell which will be extrapolated:                 eflag=-2,2  |
  |  cell which will not be extrapolated (constant):  eflag=+1,0  |
  +-------------------------------------------------------------*/

  /*-------------------------------------------+
  |  advect sca                                |
  |  time scheme: euler implicit               |
  |  advection scheme: 1st order upwind        |
  |  iteration scheme: symmetric Gauss-Seidel  |
  +-------------------------------------------*/
  const real dtau = dxmin;
  real isgn = real(iext); /* normal vector pointing to gas */
  //int mmax=8;
  int mmax=20;

  // extrapolate from eflag=1 to eflag=-2,2
  /* pseudo-time loop */
  for(int mstep=1; mstep<=mmax; mstep++){
    /* right hand side */
    for_ijk(i,j,k) {
      
      if(fabs(eflag[i][j][k])>1.0){

        real nxc = isgn*mx[i][j][k];
        real nyc = isgn*my[i][j][k];
        real nzc = isgn*mz[i][j][k];

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

#if 0
        /* add the divergence of tangential velocity */
        stmp2[i][j][k] -= rhs[i][j][k];

        /* add the effect of curvature */
        stmp2[i][j][k] -= kappa[i][j][k]*sca[i][j][k];
#endif


      }
    }
    for_aijk(i,j,k)
      stmp3[i][j][k]=0.0;

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
        if(fabs(eflag[i][j][k])>1.0){
          real diag,drhs;
          real nxc = isgn*mx[i][j][k];
          real nyc = isgn*my[i][j][k];
          real nzc = isgn*mz[i][j][k];

          diag =1.0/dtau;

#if 0
        /* add the effect of curvature */
        diag += kappa[i][j][k];
#endif

          diag+= 0.5*(nxc+fabs(nxc))/dxw(i);
          drhs = 0.5*(nxc+fabs(nxc))*stmp3[i-1][j][k]/dxw(i);

          diag-= 0.5*(nxc-fabs(nxc))/dxe(i);
          drhs-= 0.5*(nxc-fabs(nxc))*stmp3[i+1][j][k]/dxe(i);

          diag+= 0.5*(nyc+fabs(nyc))/dys(j);
          drhs+= 0.5*(nyc+fabs(nyc))*stmp3[i][j-1][k]/dys(j);

          diag-= 0.5*(nyc-fabs(nyc))/dyn(j);
          drhs-= 0.5*(nyc-fabs(nyc))*stmp3[i][j+1][k]/dyn(j);

          diag+= 0.5*(nzc+fabs(nzc))/dzb(k);
          drhs+= 0.5*(nzc+fabs(nzc))*stmp3[i][j][k-1]/dzb(k);

          diag-= 0.5*(nzc-fabs(nzc))/dzt(k);
          drhs-= 0.5*(nzc-fabs(nzc))*stmp3[i][j][k+1]/dzt(k);

          stmp3[i][j][k] = (stmp2[i][j][k]+drhs)/diag;
          //boil::oout<<i<<" "<<j<<" "<<k<<" "<<stmp2[i][j][k]<<" "<<stmp3[i][j][k]<<" "<<drhs<<" "<<diag<<boil::endl;
        }
      }}}
    }

    /* update sca */
    for_ijk(i,j,k) {
      if(fabs(eflag[i][j][k])>1.0){
        sca[i][j][k] += stmp3[i][j][k];
      }
    }
    sca.exchange();
     
    //exit(0);
    //boil::oout<<"HERE!"<<boil::endl;
  }

  return;
}

