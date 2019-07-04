#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     7x3x3 stencil
*     J.Lopez et al., Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*
*     modification: selection of cells where curvature is calculated is 
*                   modified from Lopez et al.
*
*     output: kappa
*******************************************************************************/
  /*-----------------------------------------------------+
  |  calculate curvature only where iflag=1              |
  |  Step 1: calculate normal vector                     |
  |  Step 2: calculate alpha using VOF approach          |
  |  Step 3: extrapolate alpha in wall                   |
  |  Step 4: calculate area density using marching cube  |
  |  Step 5: define iflag                                |
  +-----------------------------------------------------*/
#if 0
  for_aijk(i,j,k) {
    nx[i][j][k] = -nx[i][j][k];
    ny[i][j][k] = -ny[i][j][k];
    nz[i][j][k] = -nz[i][j][k];
    phi[i][j][k] = 1.-phi[i][j][k];
  }
#endif

  iflag=0;
#if 0  // don't use this. Yohei
  /* calculate alpha in cells */
  norm_cc(phi);     // nx, ny, nz are necessary for extract_alpha
  extract_alpha();  // this is necesarry for update_at_walls

  /* prerequisite for marching cubes */
  update_at_walls();
  //boil::plot->plot(phi, "phi", time->current_step());

  /* calculate area density */
  cal_adens();  // it access nx, ny and nz
  for_ijk(i,j,k) {
    if (adens[i][j][k]>0.0) iflag[i][j][k]=1;
  }
#else
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      if((phi[i-1][j][k]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
      if((phi[i+1][j][k]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
      if((phi[i][j-1][k]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
      if((phi[i][j+1][k]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
      if((phi[i][j][k-1]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
      if((phi[i][j][k+1]-0.5)*(phi[i][j][k]-0.5)<=0.0){ iflag[i][j][k]=1;}
    }
  }
#endif
  iflag.exchange();

#if 0
  if(boil::cart.iam()==6){
    std::cout<<"curv_HF:proc6: "<<phi[si()-1+32][sj()-1+1][sk()-1+0]<<" "
                                <<phi[si()-1+32][sj()-1+1][sk()-1+1]<<"\n";
    std::cout<<"curv_HF:iflag: "<<iflag[si()-1+32][sj()-1+1][sk()-1+0]<<" "
                                <<iflag[si()-1+32][sj()-1+1][sk()-1+1]<<"\n";
  }
#endif

  /* calculate normal vector */
  //true_norm_vect();
#if 0
  norm_young(phi);  // Young's method is good for low resolution
#else
  for_avijk(phi,i,j,k) {
    nx[i][j][k] = mx[i][j][k];
    ny[i][j][k] = my[i][j][k];
    nz[i][j][k] = mz[i][j][k];
  }
#endif

  /*------------------------+
  |  curvature calculation  |
  +------------------------*/
  //kappa=boil::unreal;
  //kappa=0.0;  // use iflag for marking cal/non-cal of body-force
  kappa=boil::unreal;


#if 0
  int ist=phi.si()-1;
  int jst=phi.sj()-1;
  int kst=phi.sk()-1;
  std::cout<<"ist="<<ist<<" "<<jst<<" "<<kst<<"\n";
  std::cout<<"iflag & kappa "<<iflag[ist+23][jst+1][kst+22]<<" "
                            <<kappa[ist+23][jst+1][kst+22]<<" "
                            <<phi[ist+23][jst+1][kst+22]<<"\n";
#endif

  /* copy phi to stmp */
  for_aijk(i,j,k) {
    stmp[i][j][k] = phi[i][j][k];
  }

  for_ijk(i,j,k) {

    /* calculate curvature only when iflag=1 */
    if (iflag[i][j][k]==1) {

      /* select direction of stencil-7 */
      int dirMax=0;
      if (fabs(nx[i][j][k])<fabs(ny[i][j][k])) {
        if (fabs(ny[i][j][k])<fabs(nz[i][j][k])) {
          dirMax=3;
        } else {
          dirMax=2;
        }
      } else {
        if (fabs(nx[i][j][k])<fabs(nz[i][j][k])) {
          dirMax=3;
        } else {
          dirMax=1;
        }
      }
      
#if 0
      if(i==5&&j==19)
        boil::oout<<"dirmax "<<i<<" "<<j<<" "<<k<<" "<<dirMax<<" | "<<nx[i][j][k]<<" "<<ny[i][j][k]<<" "<<nz[i][j][k]<<" | "<<mx[i][j][k]<<" "<<my[i][j][k]<<" "<<mz[i][j][k]<<boil::endl;
#endif

      real hmm = 0.0, hcm = 0.0, hpm = 0.0;
      real hmc = 0.0, hcc = 0.0, hpc = 0.0;
      real hmp = 0.0, hcp = 0.0, hpp = 0.0;
      real hc_limit = 0.0;

      if (dirMax==1) {
        if(!ifull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with dirMax="<<dirMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int imin=-3;
        int imax=3;  // normal stencil size 7 (=3+1+3)
        if (iminc) imin=max(-3,si()-i);  // limit stencil size for cut-stencil
        if (imaxc) imax=min( 3,ei()-i);
        if (dom->ibody().off(i-1,j,k)) { imin=-1; }  // use wall adjacent phi
        else if(dom->ibody().off(i-2,j,k)) { imin=-2; } // in solid
        if (dom->ibody().off(i+1,j,k)) { imax=1; }
        else if(dom->ibody().off(i+2,j,k)) { imax=2; }

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int ii=-1; ii>=imin; ii--) {
          if(imin==0)std::cout<<"enter\n";
          real nxc = nx[i][j][k];
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii+1][j+jj][k+kk]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nxc));
            }
          }}
          if (nxc<0.0) {
            hc_limit += phi.dxc(i+ii);
          }
        }

        for (int ii=1; ii<=imax; ii++) {
          real nxc = nx[i][j][k];
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii-1][j+jj][k+kk]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nxc));
            }
          }}
          if (nxc>0.0) {
            hc_limit += phi.dxc(i+ii);
          }
        }

        for (int ii=imin; ii<=imax; ii++) {
          hmm += stmp[i+ii][j-1][k-1]*phi.dxc(i+ii);
          hcm += stmp[i+ii][j  ][k-1]*phi.dxc(i+ii);
          hpm += stmp[i+ii][j+1][k-1]*phi.dxc(i+ii);
          hmc += stmp[i+ii][j-1][k  ]*phi.dxc(i+ii);
          hcc += stmp[i+ii][j  ][k  ]*phi.dxc(i+ii);
          hpc += stmp[i+ii][j+1][k  ]*phi.dxc(i+ii);
          hmp += stmp[i+ii][j-1][k+1]*phi.dxc(i+ii);
          hcp += stmp[i+ii][j  ][k+1]*phi.dxc(i+ii); 
          hpp += stmp[i+ii][j+1][k+1]*phi.dxc(i+ii);
        }

        //if (3*phi.dxc(i)<hcc && hcc<(4.0)*phi.dxc(i)) {
        if (hc_limit<=hcc && hcc<=(hc_limit+phi.dxc(i))) {
        real hy  = jfull*(hpc-hmc)/(dys(j)+dyn(j));
        real hz  = kfull*(hcp-hcm)/(dzb(k)+dzt(k));
        real hyy = jfull*(hpc-2.0*hcc+hmc)/(phi.dyc(j)*phi.dyc(j));
        real hzz = kfull*(hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hyz = jfull*kfull*(hpp-hpm-hmp+hmm) / (4.0*phi.dyc(j)*phi.dzc(k));
        kappa[i][j][k] = -1.0
                       * (hyy + hzz + hyy*hz*hz + hzz*hy*hy - 2.0*hyz*hy*hz)
                       / pow(1.0 + hy*hy + hz*hz, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=imin; ii<=imax; ii++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = phi[i+ii][j+jj][k+kk];
        }}}

      } else if (dirMax==2) {
        if(!jfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with dirMax="<<dirMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int jmin=-3;
        int jmax=3;  // normal stencil size 7 (=3+1+3)
        if (jminc) jmin=max(-3,sj()-j);  // limit stencil size for cut-stencil
        if (jmaxc) jmax=min( 3,ej()-j);
        //std::cout<<"j,jmin,jmax=: "<<j<<" "<<jmin<<" "<<jmax<<"\n";
        if (dom->ibody().off(i,j-1,k)) { jmin=-1; }
        else if(dom->ibody().off(i,j-2,k)) { jmin=-2; }
        if (dom->ibody().off(i,j+1,k)) { jmax=1; }
        else if(dom->ibody().off(i,j+2,k)) { jmax=2; }

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int jj=-1; jj>=jmin; jj--) {
          real nyc = ny[i][j][k];
          for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc)
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj+1][k+kk]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nyc));
            }
          }}
          if (nyc<0.0) {
            hc_limit += phi.dyc(j+jj);
          }
        }

        for (int jj=1; jj<=jmax; jj++) {
          real nyc = ny[i][j][k];
         for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc)
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj-1][k+kk]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nyc));
            }
          }}
          if (nyc>0.0) {
            hc_limit += phi.dyc(j+jj);
          }
        }

        for (int jj=jmin; jj<=jmax; jj++) {
          hmm += stmp[i-1][j+jj][k-1]*phi.dyc(j+jj);
          hcm += stmp[i  ][j+jj][k-1]*phi.dyc(j+jj);
          hpm += stmp[i+1][j+jj][k-1]*phi.dyc(j+jj);
          hmc += stmp[i-1][j+jj][k  ]*phi.dyc(j+jj);
          hcc += stmp[i  ][j+jj][k  ]*phi.dyc(j+jj);
          hpc += stmp[i+1][j+jj][k  ]*phi.dyc(j+jj);
          hmp += stmp[i-1][j+jj][k+1]*phi.dyc(j+jj);
          hcp += stmp[i  ][j+jj][k+1]*phi.dyc(j+jj);
          hpp += stmp[i+1][j+jj][k+1]*phi.dyc(j+jj);
        }

        if (hc_limit<=hcc && hcc<=(hc_limit+phi.dyc(j))) {
        //if (3*phi.dyc(j)<hcc && hcc<(4.0)*phi.dyc(j)) {
        real hx  = ifull*(hpc-hmc)/(dxw(i)+dxe(i));
        real hz  = kfull*(hcp-hcm)/(dzb(k)+dzt(k));
        real hxx = ifull*(hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hzz = kfull*(hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hxz = ifull*kfull*(hpp-hpm-hmp+hmm) / (4.0*phi.dxc(i)*phi.dzc(k));
        kappa[i][j][k] = -1.0
                       * (hxx + hzz + hxx*hz*hz + hzz*hx*hx - 2.0*hxz*hx*hz)
                       / pow(1.0 + hx*hx + hz*hz, 1.5);

        //boil::oout<<j<<" "<<hx<<" "<<hz<<" "<<hxx<<" "<<hzz<<" "<<hxz<<" "<<kappa[i][j][k]<<" | "<<stmp[i-1][j][k]<<" "<<stmp[i][j][k]<<" "<<stmp[i+1][j][k]<<boil::endl;
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=-1; ii<=1; ii++) {
        for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = phi[i+ii][j+jj][k+kk];
        }}}

      } else if (dirMax==3) {
        if(!kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with dirMax="<<dirMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int kmin=-3;
        int kmax=3;  // normal stencil size 7 (=3+1+3)
        if (kminc) kmin=max(-3,sk()-k);  // limit stencil size for cut-stencil
        if (kmaxc) kmax=min( 3,ek()-k);

        if (dom->ibody().off(i,j,k-1)) { kmin=-1; }
        else if(dom->ibody().off(i,j,k-2)) { kmin=-2; }
        if (dom->ibody().off(i,j,k+1)) { kmax=1; }
        else if(dom->ibody().off(i,j,k+2)) { kmax=2; }
#if 0
        if (boil::cart.iam()==0){
          if(i==si()-1+5 && j==sj()-1+1 && k==sk()-1+64){
            std::cout<<"0:dirMax==3 "<<kmin<<" "<<kmax<<" "<<k<<" "<<ek()<<"\n";
          }
        }
        if (boil::cart.iam()==1){
          if(i==si()-1+5 && j==sj()-1+1 && k==sk()-1+1){
            std::cout<<"1:dirMax==3 "<<kmin<<" "<<kmax<<" "<<k<<" "<<ek()<<"\n";
          }
        }
#endif

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int kk=-1; kk>=kmin; kk--) {
          real nzc = nz[i][j][k];
#if 0
        if (boil::cart.iam()==0){
          if(i==si()-1+5 && j==sj()-1+1 && k==sk()-1+64){
            std::cout<<"proc0:nzc= "<<kk<<" "<<nzc<<"\n";
          }
        }
#endif
          for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj][k+kk+1]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nzc));
            }
          }}
          if (nzc<0.0) {
            hc_limit += phi.dzc(k+kk);
          }
        }

        for (int kk=1; kk<=kmax; kk++) { 
          real nzc = nz[i][j][k];
         for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj][k+kk-1]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nzc));
            }
          }}
          if (nzc>0.0) {
            hc_limit += phi.dzc(k+kk);
          }
        }

        for (int kk=kmin; kk<=kmax; kk++) {    
          hmm += stmp[i-1][j-1][k+kk]*phi.dzc(k+kk);
          hcm += stmp[i  ][j-1][k+kk]*phi.dzc(k+kk);
          hpm += stmp[i+1][j-1][k+kk]*phi.dzc(k+kk);
          hmc += stmp[i-1][j  ][k+kk]*phi.dzc(k+kk);
          hcc += stmp[i  ][j  ][k+kk]*phi.dzc(k+kk);
          hpc += stmp[i+1][j  ][k+kk]*phi.dzc(k+kk);
          hmp += stmp[i-1][j+1][k+kk]*phi.dzc(k+kk); 
          hcp += stmp[i  ][j+1][k+kk]*phi.dzc(k+kk);
          hpp += stmp[i+1][j+1][k+kk]*phi.dzc(k+kk);
        }
#if 0
        if (boil::cart.iam()==0){
          if(i==si()-1+5 && j==sj()-1+1 && k==sk()-1+64){
            std::cout<<"check proc0:"<<kmin<<" "<<kmax<<" "<<k<<" "<<ek()<<"\n";
            std::cout<<hc_limit<<" "<<hcc<<" "<<hc_limit+phi.dzc(k)<<"\n";
          }
        }
#endif
        if (hc_limit<hcc && hcc<(hc_limit+phi.dzc(k))) {
        //if (3*phi.dzc(k)<hcc && hcc<(4.0)*phi.dzc(k)) {
        real hx  = ifull*(hpc-hmc)/(dxw(i)+dxe(i));
        real hy  = jfull*(hcp-hcm)/(dys(j)+dyn(j));
        real hxx = ifull*(hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hyy = jfull*(hcp-2.0*hcc+hcm)/(phi.dyc(j)*phi.dyc(j));
        real hxy = ifull*jfull*(hpp-hpm-hmp+hmm) / (4.0*phi.dxc(i)*phi.dyc(j));
        kappa[i][j][k] = -1.0
                       * (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.0*hxy*hx*hy)
                       / pow(1.0 + hx*hx + hy*hy, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=-1; ii<=1; ii++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {
            stmp[i+ii][j+jj][k+kk] = phi[i+ii][j+jj][k+kk];
        }}}

      }
    }
  }
  kappa.bnd_update();
  kappa.exchange();
  iflag.exchange();

#if 0
  int iproc=boil::cart.iam();
  if (iproc==0){
    std::cout<<"curv_HF: "<<kappa[si()-1+5][sj()-1+1][sk()-1+64]<<"\n";
    std::cout<<"proc0:j=0 "<<phi[si()-1+5][sj()-1+0][sk()-1+64]<<" "
                           <<phi[si()-1+5][sj()-1+0][sk()-1+65]<<" "
                           <<phi[si()-1+5][sj()-1+0][sk()-1+66]<<" "
                           <<phi[si()-1+5][sj()-1+0][sk()-1+67]<<"\n";
    std::cout<<"proc0:j=1 "<<phi[si()-1+5][sj()-1+1][sk()-1+64]<<" "
                           <<phi[si()-1+5][sj()-1+1][sk()-1+65]<<" "
                           <<phi[si()-1+5][sj()-1+1][sk()-1+66]<<" "
                           <<phi[si()-1+5][sj()-1+1][sk()-1+67]<<"\n";
    std::cout<<"proc0:j=2 "<<phi[si()-1+5][sj()-1+2][sk()-1+64]<<" "
                           <<phi[si()-1+5][sj()-1+2][sk()-1+65]<<" "
                           <<phi[si()-1+5][sj()-1+2][sk()-1+66]<<" "
                           <<phi[si()-1+5][sj()-1+2][sk()-1+67]<<"\n";
    std::cout<<"proc0:dzc "<<phi.dzc(sk()-1+64)<<" "
                           <<phi.dzc(sk()-1+65)<<" "
                           <<phi.dzc(sk()-1+66)<<" "
                           <<phi.dzc(sk()-1+67)<<"\n";
  }
  if (iproc==1){
    std::cout<<"curv_HF: "<<kappa[si()-1+5][sj()-1+1][sk()-1+1]<<"\n";
  }
  exit(0);
#endif

#if 0
  //if(time->current_step()==2) {
  /* visualize iflag */
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=iflag[i][j][k];
  }
  boil::plot->plot(phi,kappa,stmp, "clr-kappa-iflag", time->current_step());
  exit(0);
  //} 
#endif
#if 0
  int iproc=boil::cart.iam();
  if (iproc==0){
    std::cout<<"curv_HF: "<<kappa[si()-1+5][sj()-1+1][sk()-1+64]<<"\n";
  }
  if (iproc==1){
    std::cout<<"curv_HF: "<<kappa[si()-1+5][sj()-1+1][sk()-1+1]<<"\n";
  }
#endif

#if 1
  // extrapolate kappa
  stmp = kappa;
  iflagx = iflag;

  //for(int iloop=1; iloop<2; iloop++) {
  for(int iloop=1; iloop<3; iloop++) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(iflag[i][j][k]==0) {
        int inb =  iflag[i-1][j][k] + iflag[i+1][j][k]
                 + iflag[i][j-1][k] + iflag[i][j+1][k]
                 + iflag[i][j][k-1] + iflag[i][j][k+1];
        if (inb >= 1) {
            stmp[i][j][k] = (real(iflag[i-1][j][k]) * kappa[i-1][j][k]
                           + real(iflag[i+1][j][k]) * kappa[i+1][j][k]
                           + real(iflag[i][j-1][k]) * kappa[i][j-1][k]
                           + real(iflag[i][j+1][k]) * kappa[i][j+1][k]
                           + real(iflag[i][j][k-1]) * kappa[i][j][k-1]
                           + real(iflag[i][j][k+1]) * kappa[i][j][k+1])
                           /real(inb);
            iflagx[i][j][k] = 2;  // iflag=2 for extrapolated
        }
      }
    }
    stmp.exchange();
    iflagx.exchange();
    kappa = stmp;
    iflag = iflagx;
  }
#endif

  bdcurv();

#if 0
  for_ijk(i,j,k) {
    if(iflag[i][j][k]==1&&fabs(kappa[i][j][k])<700.)
      boil::oout<<"test "<<i<<" "<<j<<" "<<k<<" "<<kappa[i][j][k]<<boil::endl;
  }
  //if(time->current_step()==1) {
  //if(time->current_step()%100==0) {
  /* visualize iflag */
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=iflag[i][j][k];
  }
  boil::plot->plot(phi,kappa,stmp, "clr-kappa-iflag", time->current_step());
  //exit(0);
  //} 
#endif

#if 0
  for_aijk(i,j,k) {
    nx[i][j][k] = -nx[i][j][k];
    ny[i][j][k] = -ny[i][j][k];
    nz[i][j][k] = -nz[i][j][k];
    phi[i][j][k] = 1.-phi[i][j][k];
    //if(boil::realistic(kappa[i][j][k]))
    //    kappa[i][j][k] = -kappa[i][j][k];
  }
#endif

  return;
}

