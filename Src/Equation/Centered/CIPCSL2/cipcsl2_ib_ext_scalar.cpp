#include "cipcsl2.h"
#ifdef IB
//#define DEBUG
using namespace std;

/******************************************************************************/
void CIPCSL2::ib_ext_scalar(Scalar & sca) {
/***************************************************************************//**
*  \brief Set immersed boundary value of color function
*******************************************************************************/
  if(dom->ibody().ncall()==0) return;

#if 0
  /* set iflag for immersed boundary */
  ib_set_iflag();
#ifdef DEBUG
  std::cout<<"ib_ext_color:ib_st_iflag \n";
#endif


  for(int loop=0; loop<16; loop++){  
    for_ijk(i,j,k){

      // calculate only solid next to fluid
      if(iflag[i][j][k]==-1){

        // i-direction
        real uc, dsdxm, dsdxp;
        uc = dom->ibody().nwx(i,j,k);
        dsdxm = (sca[i][j][k]-sca[i-1][j][k])/sca.dxw(i);
        dsdxp = (sca[i+1][j][k]-sca[i][j][k])/sca.dxe(i);

        // j-direction
        real vc, dsdym, dsdyp;
        vc = dom->ibody().nwy(i,j,k);
        dsdym = (sca[i][j][k]-sca[i][j-1][k])/sca.dys(j);
        dsdyp = (sca[i][j+1][k]-sca[i][j][k])/sca.dyn(j);

        // k-direction
        real wc, dsdzm, dsdzp;
        wc = dom->ibody().nwz(i,j,k);
        dsdzm = (sca[i][j][k]-sca[i][j][k-1])/sca.dzb(k);
        dsdzp = (sca[i][j][k+1]-sca[i][j][k])/sca.dzt(k);

        fn[i][j][k] = 0.5*(uc+fabs(uc))*dsdxm
                        + 0.5*(uc-fabs(uc))*dsdxp
                        + 0.5*(vc+fabs(vc))*dsdym
                        + 0.5*(vc-fabs(vc))*dsdyp
                        + 0.5*(wc+fabs(wc))*dsdzm
                        + 0.5*(wc-fabs(wc))*dsdzp;
      }
    }

    // update
    for_ijk(i,j,k){
      // update only solid next to fluid
      if(iflag[i][j][k]==-1){
        int iof,jof,kof;
        iof=jof=kof=0;
        if(dom->ibody().nwz(i,j,k)<-0.7)     {kof= 1;}  // crude code
        else if(dom->ibody().nwz(i,j,k)> 0.7) {kof=-1;}
        else if(dom->ibody().nwx(i,j,k)<-0.7) {iof= 1;}
        else if(dom->ibody().nwx(i,j,k)< 0.7) {iof= 1;}
        else if(dom->ibody().nwy(i,j,k)<-0.7) {jof= 1;}
        else if(dom->ibody().nwy(i,j,k)< 0.7) {jof= 1;}
        real dmin=boil::minr(sca.dxc(i+iof),sca.dyc(j+jof),sca.dzc(k+kof));
	real cfl = 0.5;
        real dtau = cfl*dmin;
        sca[i][j][k] -= dtau * fn[i][j][k];
      }
    }
    sca.bnd_update();
    sca.exchange();
  }

#if 0
    Scalar ux(*phi.domain());
    Scalar uy(*phi.domain());
    Scalar uz(*phi.domain());
    for_ijk(i,j,k) {
      ux[i][j][k]=dom->ibody().nwx(i,j,k);
      uy[i][j][k]=dom->ibody().nwy(i,j,k);
      uz[i][j][k]=dom->ibody().nwz(i,j,k);
    }
    boil::plot->plot(sca, ux, uy, uz, 
          "ib_init-sca-ux-uy-uz", time->current_step());
    exit(0);
#endif
#else

  // extrapolate clr  //crude 2014.01.16, Ver1.1.10.23R2
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=dom->ibody().nwx(i,j,k);
    real uy=dom->ibody().nwy(i,j,k);
    real uz=dom->ibody().nwz(i,j,k);
    Dir d = Dir::undefined();
    if (fabs(uz)>0.707) {
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
      }
    } else if (fabs(ux)>0.707) {
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
      }
    } else if (fabs(uy)>0.707) {
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
      }
    } else {
      std::cout<<"cipcsl2_ib_ext_scalar: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*----------------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in computational domain  |
    |         (i+iof, j+jof, k+kof) is on wall                  |
    +----------------------------------------------------------*/
    if(dom->ibody().off(i+iof,j+jof,k+kof)){
      sca[i+iof][j+jof][k+kof] = sca[i][j][k];
    }
  }
  sca.bnd_update();
  sca.exchange();
#endif

}
#endif
