#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curv_HF() {
/******************************************************************************/
  curv_HF(phi,true);
   
  return;
}

/******************************************************************************/
void VOF::curv_HF(const Scalar & Phi, const bool anc_flag) {
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

  if(anc_flag) {

  /* calculate normal vector */
  //norm_cc(Phi);
  norm_young(Phi);  // Young's method is good for low resolution

  /* calculate alpha in cells */
  extract_alpha();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate area density */
  cal_adens();

  /* normal vector */
  true_norm_vect();

  }

  iflag=0;
  for_ijk(i,j,k) {
    if (adens[i][j][k]>0.0) iflag[i][j][k]=1;
  }
  iflag.exchange();

  /*------------------------+
  |  curvature calculation  |
  +------------------------*/
  kappa=boil::unreal;
#if 0
  int ist=Phi.si()-1;
  int jst=Phi.sj()-1;
  int kst=Phi.sk()-1;
  std::cout<<"ist="<<ist<<" "<<jst<<" "<<kst<<"\n";
  std::cout<<"iflag & kappa "<<iflag[ist+23][jst+1][kst+22]<<" "
                            <<kappa[ist+23][jst+1][kst+22]<<" "
                            <<Phi[ist+23][jst+1][kst+22]<<"\n";
#endif

  /* copy Phi to stmp */
  for_aijk(i,j,k) {
    stmp[i][j][k] = Phi[i][j][k];
  }

  for_ijk(i,j,k) {

    /* calculate curvature only when iflag=1 */
    if (iflag[i][j][k]==1) {

      /* select direction of stencil-7 */
      int dirMax=0;
      if (fabs(mx[i][j][k])<fabs(my[i][j][k])) {
        if (fabs(my[i][j][k])<fabs(mz[i][j][k])) {
          dirMax=3;
        } else {
          dirMax=2;
        }
      } else {
        if (fabs(mx[i][j][k])<fabs(mz[i][j][k])) {
          dirMax=3;
        } else {
          dirMax=1;
        }
      }

      real hmm = 0.0, hcm = 0.0, hpm = 0.0;
      real hmc = 0.0, hcc = 0.0, hpc = 0.0;
      real hmp = 0.0, hcp = 0.0, hpp = 0.0;
      real hc_limit = 0.0;

      if (dirMax==1) {
        /* check stencil size */
#if 0
        int imin=max(-3,si()-i);
        int imax=min( 3,ei()-i);
#else
        int imin(-3), imax(3);
#endif

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int ii=-1; ii>=imin; ii--) {
          if(imin==0)std::cout<<"enter\n";
          real nxc = mx[i][j][k];
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii+1][j+jj][k+kk]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nxc));
            }
          }}
          if (nxc<0.0) {
            hc_limit += Phi.dxc(i+ii);
          }
        }

        for (int ii=1; ii<=imax; ii++) {
          real nxc = mx[i][j][k];
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii-1][j+jj][k+kk]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nxc));
            }
          }}
          if (nxc>0.0) {
            hc_limit += Phi.dxc(i+ii);
          }
        }

        for (int ii=imin; ii<=imax; ii++) {
          hmm += stmp[i+ii][j-1][k-1]*Phi.dxc(i+ii);
          hcm += stmp[i+ii][j  ][k-1]*Phi.dxc(i+ii);
          hpm += stmp[i+ii][j+1][k-1]*Phi.dxc(i+ii);
          hmc += stmp[i+ii][j-1][k  ]*Phi.dxc(i+ii);
          hcc += stmp[i+ii][j  ][k  ]*Phi.dxc(i+ii);
          hpc += stmp[i+ii][j+1][k  ]*Phi.dxc(i+ii);
          hmp += stmp[i+ii][j-1][k+1]*Phi.dxc(i+ii);
          hcp += stmp[i+ii][j  ][k+1]*Phi.dxc(i+ii); 
          hpp += stmp[i+ii][j+1][k+1]*Phi.dxc(i+ii);
        }

        //if (3*Phi.dxc(i)<hcc && hcc<(4.0)*Phi.dxc(i)) {
        if (hc_limit<=hcc && hcc<=(hc_limit+Phi.dxc(i))) {
        real hy  = (hpc-hmc)/(dys(j)+dyn(j));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hyy = (hpc-2.0*hcc+hmc)/(Phi.dyc(j)*Phi.dyc(j));
        real hzz = (hcp-2.0*hcc+hcm)/(Phi.dzc(k)*Phi.dzc(k));
        real hyz = (hpp-hpm-hmp+hmm) / (4.0*Phi.dyc(j)*Phi.dzc(k));
        kappa[i][j][k] = -1.0
                       * (hyy + hzz + hyy*hz*hz + hzz*hy*hy - 2.0*hyz*hy*hz)
                       / pow(1.0 + hy*hy + hz*hz, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=imin; ii<=imax; ii++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = Phi[i+ii][j+jj][k+kk];
        }}}

      } else if (dirMax==2) {
        /* check stencil size */
#if 0
        int jmin=max(-3,sj()-j);
        int jmax=min( 3,ej()-j);
#else
        int jmin(-3), jmax(3);
#endif
        //std::cout<<"j,jmin,jmax=: "<<j<<" "<<jmin<<" "<<jmax<<"\n";

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int jj=-1; jj>=jmin; jj--) {
          real nyc = my[i][j][k];
          for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc)
                *(Phi[i+ii][j+jj][k+kk]-Phi[i+ii][j+jj+1][k+kk]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nyc));
            }
          }}
          if (nyc<0.0) {
            hc_limit += Phi.dyc(j+jj);
          }
        }

        for (int jj=1; jj<=jmax; jj++) {
          real nyc = my[i][j][k];
         for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc)
                *(Phi[i+ii][j+jj][k+kk]-Phi[i+ii][j+jj-1][k+kk]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nyc));
            }
          }}
          if (nyc>0.0) {
            hc_limit += Phi.dyc(j+jj);
          }
        }

        for (int jj=jmin; jj<=jmax; jj++) {
          hmm += stmp[i-1][j+jj][k-1]*Phi.dyc(j+jj);
          hcm += stmp[i  ][j+jj][k-1]*Phi.dyc(j+jj);
          hpm += stmp[i+1][j+jj][k-1]*Phi.dyc(j+jj);
          hmc += stmp[i-1][j+jj][k  ]*Phi.dyc(j+jj);
          hcc += stmp[i  ][j+jj][k  ]*Phi.dyc(j+jj);
          hpc += stmp[i+1][j+jj][k  ]*Phi.dyc(j+jj);
          hmp += stmp[i-1][j+jj][k+1]*Phi.dyc(j+jj);
          hcp += stmp[i  ][j+jj][k+1]*Phi.dyc(j+jj);
          hpp += stmp[i+1][j+jj][k+1]*Phi.dyc(j+jj);
        }

        if (hc_limit<=hcc && hcc<=(hc_limit+Phi.dyc(j))) {
        //if (3*Phi.dyc(j)<hcc && hcc<(4.0)*Phi.dyc(j)) {
        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hxx = (hpc-2.0*hcc+hmc)/(Phi.dxc(i)*Phi.dxc(i));
        real hzz = (hcp-2.0*hcc+hcm)/(Phi.dzc(k)*Phi.dzc(k));
        real hxz = (hpp-hpm-hmp+hmm) / (4.0*Phi.dxc(i)*Phi.dzc(k));
        kappa[i][j][k] = -1.0
                       * (hxx + hzz + hxx*hz*hz + hzz*hx*hx - 2.0*hxz*hx*hz)
                       / pow(1.0 + hx*hx + hz*hz, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=-1; ii<=1; ii++) {
        for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = Phi[i+ii][j+jj][k+kk];
        }}}

      } else if (dirMax==3) {
        /* check stencil size */
#if 0
        int kmin=max(-3,sk()-k);
        int kmax=min( 3,ek()-k);
#else
        int kmin(-3), kmax(3);
#endif
        //std::cout<<"k,kmin,kmax=: "<<k<<" "<<kmin<<" "<<kmax<<"\n";

        /* Calculate Eq. (2) in J.Lopez et al. */
        for (int kk=-1; kk>=kmin; kk--) {
          real nzc = mz[i][j][k];
          for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)
                *(Phi[i+ii][j+jj][k+kk]-Phi[i+ii][j+jj][k+kk+1]) > 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nzc));
            }
          }}
          if (nzc<0.0) {
            hc_limit += Phi.dzc(k+kk);
          }
        }

        for (int kk=1; kk<=kmax; kk++) { 
          real nzc = mz[i][j][k];
         for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)
                *(Phi[i+ii][j+jj][k+kk]-Phi[i+ii][j+jj][k+kk-1]) < 0.0){
              stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nzc));
            }
          }}
          if (nzc>0.0) {
            hc_limit += Phi.dzc(k+kk);
          }
        }

        for (int kk=kmin; kk<=kmax; kk++) {    
          hmm += stmp[i-1][j-1][k+kk]*Phi.dzc(k+kk);
          hcm += stmp[i  ][j-1][k+kk]*Phi.dzc(k+kk);
          hpm += stmp[i+1][j-1][k+kk]*Phi.dzc(k+kk);
          hmc += stmp[i-1][j  ][k+kk]*Phi.dzc(k+kk);
          hcc += stmp[i  ][j  ][k+kk]*Phi.dzc(k+kk);
          hpc += stmp[i+1][j  ][k+kk]*Phi.dzc(k+kk);
          hmp += stmp[i-1][j+1][k+kk]*Phi.dzc(k+kk); 
          hcp += stmp[i  ][j+1][k+kk]*Phi.dzc(k+kk);
          hpp += stmp[i+1][j+1][k+kk]*Phi.dzc(k+kk);
        }

        if (hc_limit<=hcc && hcc<=(hc_limit+Phi.dzc(k))) {
        //if (3*Phi.dzc(k)<hcc && hcc<(4.0)*Phi.dzc(k)) {
        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hy  = (hcp-hcm)/(dys(j)+dyn(j));
        real hxx = (hpc-2.0*hcc+hmc)/(Phi.dxc(i)*Phi.dxc(i));
        real hyy = (hcp-2.0*hcc+hcm)/(Phi.dyc(j)*Phi.dyc(j));
        real hxy = (hpp-hpm-hmp+hmm) / (4.0*Phi.dxc(i)*Phi.dyc(j));
        kappa[i][j][k] = -1.0
                       * (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.0*hxy*hx*hy)
                       / pow(1.0 + hx*hx + hy*hy, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=-1; ii<=1; ii++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int kk=kmin; kk<=kmax; kk++) {
            stmp[i+ii][j+jj][k+kk] = Phi[i+ii][j+jj][k+kk];
        }}}

      }
    }
  }

#if 0
  for_vijk(kappa,i,j,k) {
    if(iflag[i][j][k]!=0) {
      kappa[i][j][k] = 1.001/1e-3;
    }
  }
  diffuse(kappa,stmp4,4);
  for_avijk(kappa,i,j,k) {
    //stmp2[i][j][k] = kappa[i][j][k];
    kappa[i][j][k] = stmp4[i][j][k];
  }
#endif

  kappa.exchange();
  iflag.exchange();

#if 0
  //if(time->current_step()==2) {
  if(!anc_flag) {
  /* visualize iflag */
  boil::plot->plot(Phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=iflag[i][j][k];
  }
  boil::plot->plot(Phi,kappa,stmp, "clr-kappa-iflag", time->current_step());
  exit(0);
  } 
#endif

#if 1
  // extrapolate kappa
  stmp = kappa;
  iflagx = iflag;

  //for(int iloop=1; iloop<2; iloop++) {
  for(int iloop=1; iloop<3; iloop++) {
    for_ijk(i,j,k) {
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
            iflagx[i][j][k] = 1;
        }
      }
    }
    stmp.exchange();
    iflagx.exchange();
    kappa = stmp;
    iflag = iflagx;
  }
#endif

#if 0
  if(time->current_step()==1) {
  /* visualize iflag */
  boil::plot->plot(Phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=iflag[i][j][k];
  }
  boil::plot->plot(Phi,kappa,stmp, "clr-kappa-iflag", time->current_step());
  exit(0);
  } 
#endif
  return;
}

