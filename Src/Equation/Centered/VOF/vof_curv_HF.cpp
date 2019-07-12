#include "vof.h"
#define LOPEZ
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
  const real theta_crit = 0.8;
  /*-----------------------------------------------------------+
  |  Step 1: calculate normal vector                           |
  |  Step 2: define iflag=1                                    |
  |  Step 3: calculate curvature if (iflag==1 & 3<=height<=4)  |
  |  Step 4: extrapolate curvatuer to two layers               |
  |  Step 5: calculate curvature at wall                       |
  +-----------------------------------------------------------*/

  /* calculate normal vector */
  // Normal vector is used for (i) dirMax and (ii) Eq. (2) in J.Lopez et al.
#if 0
  norm_young(phi);  // Young's method is good for low resolution
#else
  for_aijk(i,j,k) {
    nx[i][j][k] = mx[i][j][k];
    ny[i][j][k] = my[i][j][k];
    nz[i][j][k] = mz[i][j][k];
  }
#endif

  /* define iflag */
  iflag=0;
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
  iflag.exchange();
  //boil::plot->plot(phi,nx,ny,nz, "c-nx-ny-nz", time->current_step());

  /*------------------------+
  |  curvature calculation  |
  +------------------------*/
  kappa=boil::unreal;

  /* copy phi to stmp */
  for_aijk(i,j,k) {
    stmp[i][j][k] = min(1.0,max(0.0,phi[i][j][k]));
  }

  for_ijk(i,j,k) {

    /* calculate curvature only when iflag=1 */
    if (iflag[i][j][k]==1) {

      /* select direction of stencil-7 */
      int dirMax=0;
      real abs_nx = fabs(nx[i][j][k]);
      real abs_ny = fabs(ny[i][j][k]);
      real abs_nz = fabs(nz[i][j][k]);
      real max_abs_n;
 
      if (abs_nx<abs_ny) {
        if (abs_ny<abs_nz) {
          dirMax=3;
          max_abs_n = abs_nz;
        } else {
          dirMax=2;
          max_abs_n = abs_ny;
        }
      } else {
        if (abs_nx<abs_nz) {
          dirMax=3;
          max_abs_n = abs_nz;
        } else {
          dirMax=1;
          max_abs_n = abs_nx;
        }
      }

      real hmm = 0.0, hcm = 0.0, hpm = 0.0;
      real hmc = 0.0, hcc = 0.0, hpc = 0.0;
      real hmp = 0.0, hcp = 0.0, hpp = 0.0;
      real hc_limit = 0.0;
      real hc = 0.0;

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
        // Warning! if nxc is low accuracy (wrong sign), this algorithm
        // doesn't work !!!  (Yohei)
#ifdef LOPEZ
        real nxc = nx[i][j][k];
        for (int ii=-1; ii>=imin; ii--) {
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)  // Lopez
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii+1][j+jj][k+kk]) > 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii+1][j+jj][k+kk]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i     ][j+jj][k+kk]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii+1][j+jj][k+kk]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nxc));  // original
              }
            }
          }}
        }

        for (int ii=1; ii<=imax; ii++) {
          for (int jj=-1; jj<=1; jj++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nxc)  // Lopez
                *(stmp[i+ii][j+jj][k+kk]-stmp[i+ii-1][j+jj][k+kk]) < 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii-1][j+jj][k+kk]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i     ][j+jj][k+kk]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii-1][j+jj][k+kk]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nxc)); // original
              }
            }
          }}
        }
#endif
        // calculate hc_limit
        if (nx[i][j][k]<0.0) {
          for (int ii=-1; ii>=imin; ii--) {
            hc_limit += 1.0;
          }
        } else {
          for (int ii=1; ii<=imax; ii++) {
            hc_limit += 1.0;
          }
        }

        // calculate height
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
          hc  += stmp[i+ii][j  ][k  ];
        }

        if (hc_limit<=hc && hc<=(hc_limit+1.0)) {
          real theta = acos(max_abs_n);
          real g = 0.0;
          if (theta < theta_crit) {   // Eq.6 Lopez (2009) Comp Methods Appl
            g = 0.2;
          }
#if 0
          real hy  = jfull*(hpc-hmc)/(dys(j)+dyn(j));
          real hz  = kfull*(hcp-hcm)/(dzb(k)+dzt(k));
          real hyy = jfull*(hpc-2.0*hcc+hmc)/(phi.dyc(j)*phi.dyc(j));
          real hzz = kfull*(hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
#else
          real hy  = jfull*(g*(hpp-hmp)+(hpc-hmc)+g*(hpm-hmm))
                    / ((dys(j)+dyn(j))*(1.0+2.0*g));
          real hz  = kfull*(g*(hpp-hpm)+(hcp-hcm)+g*(hmp-hmm))
                    / ((dzb(k)+dzt(k))*(1.0+2.0*g));
          real hyy = jfull*(g*(hpp-2.0*hcp+hmp)+(hpc-2.0*hcc+hmc)+g*(hpm-2.0*hcm+hmm))
                    / ((phi.dyc(j)*phi.dyc(j))*(1.0+2.0*g));
          real hzz = kfull*(g*(hpp-2.0*hpc+hpm)+(hcp-2.0*hcc+hcm)+g*(hmp-2.0*hmc+hmm))
                    / ((phi.dzc(k)*phi.dzc(k))*(1.0+2.0*g));
#endif
          real hyz = jfull*kfull*(hpp-hpm-hmp+hmm)/(4.0*phi.dyc(j)*phi.dzc(k));
          kappa[i][j][k] = -1.0
                         * (hyy + hzz + hyy*hz*hz + hzz*hy*hy - 2.0*hyz*hy*hz)
                         / pow(1.0 + hy*hy + hz*hz, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=imin; ii<=imax; ii++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = min(1.0,max(0.0,phi[i+ii][j+jj][k+kk]));
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
        if (dom->ibody().off(i,j-1,k)) { jmin=-1; }
        else if(dom->ibody().off(i,j-2,k)) { jmin=-2; }
        if (dom->ibody().off(i,j+1,k)) { jmax=1; }
        else if(dom->ibody().off(i,j+2,k)) { jmax=2; }

        /* Calculate Eq. (2) in J.Lopez et al. */
#ifdef LOPEZ
        real nyc = ny[i][j][k];
        for (int jj=-1; jj>=jmin; jj--) {
          for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc) // Lopez
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj+1][k+kk]) > 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii][j+jj+1][k+kk]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i   ][j+jj  ][k+kk]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii][j+jj+1][k+kk]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nyc));  // original
              }
            }
          }}
        }

        for (int jj=1; jj<=jmax; jj++) {
          for (int ii=-1; ii<=1; ii++) {
          for (int kk=-1; kk<=1; kk++) {
            if (copysign(1.0, nyc) // Lopez
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj-1][k+kk]) < 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii][j+jj-1][k+kk]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i   ][j+jj  ][k+kk]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii][j+jj-1][k+kk]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nyc));
              }
            }
          }
          }
        }
#endif
        // calculate hc_limit
        if (ny[i][j][k]<0.0) {
          for (int jj=-1; jj>=jmin; jj--) {
            hc_limit += 1.0;
          }
        } else {
          for (int jj=1; jj<=jmax; jj++) {
            hc_limit += 1.0;
          }
        }

        // calculate height
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
          hc  += stmp[i  ][j+jj][k  ];
        }

        if (hc_limit<=hc && hc<=(hc_limit+1.0)) {
          real theta = acos(max_abs_n);
          real g = 0.0;
          if (theta < theta_crit) {   // Eq.6 Lopez (2009) Comp Methods Appl
            g = 0.2;
          }
#if 0
          real hx  = ifull*(hpc-hmc)/(dxw(i)+dxe(i));
          real hz  = kfull*(hcp-hcm)/(dzb(k)+dzt(k));
          real hxx = ifull*(hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
          real hzz = kfull*(hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
#else
          real hx  = ifull*(g*(hpp-hmp)+(hpc-hmc)+g*(hpm-hmm))
                    / ((dxw(i)+dxe(i))*(1.0+2.0*g));
          real hz  = kfull*(g*(hpp-hpm)+(hcp-hcm)+g*(hmp-hmm))
                    / ((dzb(k)+dzt(k))*(1.0+2.0*g));
          real hxx = ifull*(g*(hpp-2.0*hcp+hmp)+(hpc-2.0*hcc+hmc)+g*(hpm-2.0*hcm+hmm))
                    / ((phi.dxc(i)*phi.dxc(i))*(1.0+2.0*g));
          real hzz = kfull*(g*(hpp-2.0*hpc+hpm)+(hcp-2.0*hcc+hcm)+g*(hmp-2.0*hmc+hmm))
                    / ((phi.dzc(k)*phi.dzc(k))*(1.0+2.0*g));
#endif
          real hxz = ifull*kfull*(hpp-hpm-hmp+hmm) / (4.0*phi.dxc(i)*phi.dzc(k));
          kappa[i][j][k] = -1.0
                         * (hxx + hzz + hxx*hz*hz + hzz*hx*hx - 2.0*hxz*hx*hz)
                         / pow(1.0 + hx*hx + hz*hz, 1.5);
        } else {
          iflag[i][j][k]=0;
        }

        for (int ii=-1; ii<=1; ii++) {
        for (int jj=jmin; jj<=jmax; jj++) {
        for (int kk=-1; kk<=1; kk++) {
          stmp[i+ii][j+jj][k+kk] = min(1.0,max(0.0,phi[i+ii][j+jj][k+kk]));
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

        /* Calculate Eq. (2) in J.Lopez et al. */
#ifdef LOPEZ
        real nzc = nz[i][j][k];
        for (int kk=-1; kk>=kmin; kk--) {
          for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)  // Lopez
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj][k+kk+1]) > 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii][j+jj][k+kk+1]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i   ][j+jj][k+kk  ]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii][j+jj][k+kk+1]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0-copysign(1.0, nzc)); // original
              }
            }
          }}
        }

        for (int kk=1; kk<=kmax; kk++) { 
         for (int ii=-1; ii<=1; ii++) {
          for (int jj=-1; jj<=1; jj++) {
            if (copysign(1.0, nzc)  // Lopez
                *(phi[i+ii][j+jj][k+kk]-phi[i+ii][j+jj][k+kk-1]) < 0.0){
              if((stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i+ii][j+jj][k+kk-1]-0.5)>0.0 &&
                 (stmp[i+ii][j+jj][k+kk]-0.5)*(stmp[i   ][j+jj][k+kk  ]-0.5)>0.0) {
                stmp[i+ii][j+jj][k+kk]=stmp[i+ii][j+jj][k+kk-1]; // Yohei
              } else {
                stmp[i+ii][j+jj][k+kk]=0.5*(1.0+copysign(1.0, nzc)); // original
              }
            }
          }}
        }
#endif
        // calculate hc_limit
        if (nz[i][j][k]<0.0) {
          for (int kk=-1; kk>=kmin; kk--) {
            hc_limit += 1.0;
          }
        } else {
          for (int kk=1; kk<=kmax; kk++) {
            hc_limit += 1.0;
          }
        }

        // calculate height
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
          hc  += stmp[i  ][j  ][k+kk];
        }

        if (hc_limit<hc && hc<(hc_limit+1.0)) {
          real theta = acos(max_abs_n);
          real g = 0.0;
          if (theta < theta_crit) {   // Eq.6 Lopez (2009) Comp Methods Appl
            g = 0.2;
          }
#if 0
          real hx  = ifull*(hpc-hmc)/(dxw(i)+dxe(i));
          real hy  = jfull*(hcp-hcm)/(dys(j)+dyn(j));
          real hxx = ifull*(hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
          real hyy = jfull*(hcp-2.0*hcc+hcm)/(phi.dyc(j)*phi.dyc(j));
#else
          real hx  = ifull*(g*(hpp-hmp)+(hpc-hmc)+g*(hpm-hmm))
                    / ((dxw(i)+dxe(i))*(1.0+2.0*g));
          real hy  = jfull*(g*(hpp-hpm)+(hcp-hcm)+g*(hmp-hmm))
                    / ((dys(j)+dyn(j))*(1.0+2.0*g));
          real hxx = ifull*(g*(hpp-2.0*hcp+hmp)+(hpc-2.0*hcc+hmc)+g*(hpm-2.0*hcm+hmm))
                    / ((phi.dxc(i)*phi.dxc(i))*(1.0+2.0*g));
          real hyy = jfull*(g*(hpp-2.0*hpc+hpm)+(hcp-2.0*hcc+hcm)+g*(hmp-2.0*hmc+hmm))
                    / ((phi.dyc(j)*phi.dyc(j))*(1.0+2.0*g));
#endif
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
          stmp[i+ii][j+jj][k+kk] = min(1.0,max(0.0,phi[i+ii][j+jj][k+kk]));
        }}}

      }
    }
  }
  kappa.exchange();
  iflag.exchange();

#if 1
  // extrapolate kappa
  stmp = kappa;
  iflagx = iflag;

  //for(int iloop=1; iloop<2; iloop++) {
  //for(int iloop=1; iloop<3; iloop++) {
  for(int iloop=1; iloop<5; iloop++) { // 2019.07.09
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(iflag[i][j][k]==0) {
        int inb =  min(1,iflag[i-1][j][k]) + min(1,iflag[i+1][j][k])
                 + min(1,iflag[i][j-1][k]) + min(1,iflag[i][j+1][k])
                 + min(1,iflag[i][j][k-1]) + min(1,iflag[i][j][k+1]);
        if (inb >= 1) {
            stmp[i][j][k] = (real(min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                           + real(min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                           + real(min(1,iflag[i][j-1][k])) * kappa[i][j-1][k]
                           + real(min(1,iflag[i][j+1][k])) * kappa[i][j+1][k]
                           + real(min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                           + real(min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
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
  if(time->current_step()==1) {
  //if(time->current_step()%100==0) {
  /* visualize iflag */
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=iflag[i][j][k];
  }
  boil::plot->plot(phi,kappa,stmp, "clr-kappa-iflag", time->current_step());
  exit(0);
  } 
#endif
  return;
}

