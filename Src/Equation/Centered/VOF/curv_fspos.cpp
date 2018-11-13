#include "vof.h"

/******************************************************************************/
void VOF::curv_intpos() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  /*----------------------------------+
  |  calculate free surface position  |
  +----------------------------------*/
  cal_fs();

#if 0
  /* interfacial cells: iflag=1 */
  iflag=0;
  for_ijk(i,j,k) {
    real sum_phi=0.0;
    for(int ii=-1; ii<=1; ii++) {
      for(int jj=-1; jj<=1; jj++) {
        for(int kk=-1; kk<=1; kk++) {
          sum_phi += phi[i+ii][j+jj][k+kk];
        }
      }
    }
    if ((boil::micro<sum_phi)&&(sum_phi<27.0-boil::micro)) {
      iflag[i][j][k]=1;
    }
  }
  iflag.exchange();

  // crude code: need boundary treatment for iflag //

  /* curvature calculation */
  kappa=0.0;
  for_ijk(i,j,k) {
    if (iflag[i][j][k]==1) {
      real nx = (phi[i+1][j][k]-phi[i-1][j][k])/(dxw(i)+dxe(i));
      real ny = (phi[i][j+1][k]-phi[i][j-1][k])/(dys(j)+dyn(j));
      real nz = (phi[i][j][k+1]-phi[i][j][k-1])/(dzb(k)+dzt(k));

      int dirMax=0;
      if (fabs(nx)<fabs(ny)) {
        if (fabs(ny)<fabs(nz)) {
          dirMax=3;
        } else {
          dirMax=2;
        }
      } else {
        if (fabs(nx)<fabs(nz)) {
          dirMax=3;
        } else {
          dirMax=1;
        }
      }

      
      if (dirMax==1) {
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
        hmm = phi[i-1][j-1][k-1]*phi.dxc(i-1)
            + phi[i  ][j-1][k-1]*phi.dxc(i)
            + phi[i+1][j-1][k-1]*phi.dxc(i+1);
        hcm = phi[i-1][j  ][k-1]*phi.dxc(i-1)
            + phi[i  ][j  ][k-1]*phi.dxc(i)
            + phi[i+1][j  ][k-1]*phi.dxc(i+1);
        hpm = phi[i-1][j+1][k-1]*phi.dxc(i-1)
            + phi[i  ][j+1][k-1]*phi.dxc(i)
            + phi[i+1][j+1][k-1]*phi.dxc(i+1);
        hmc = phi[i-1][j-1][k  ]*phi.dxc(i-1)
            + phi[i  ][j-1][k  ]*phi.dxc(i)
            + phi[i+1][j-1][k  ]*phi.dxc(i+1);
        hcc = phi[i-1][j  ][k  ]*phi.dxc(i-1)
            + phi[i  ][j  ][k  ]*phi.dxc(i)
            + phi[i+1][j  ][k  ]*phi.dxc(i+1);
        hpc = phi[i-1][j+1][k  ]*phi.dxc(i-1)
            + phi[i  ][j+1][k  ]*phi.dxc(i)
            + phi[i+1][j+1][k  ]*phi.dxc(i+1);
        hmp = phi[i-1][j-1][k+1]*phi.dxc(i-1)
            + phi[i  ][j-1][k+1]*phi.dxc(i)
            + phi[i+1][j-1][k+1]*phi.dxc(i+1);
        hcp = phi[i-1][j  ][k+1]*phi.dxc(i-1)
            + phi[i  ][j  ][k+1]*phi.dxc(i)
            + phi[i+1][j  ][k+1]*phi.dxc(i+1);
        hpp = phi[i-1][j+1][k+1]*phi.dxc(i-1)
            + phi[i  ][j+1][k+1]*phi.dxc(i)
            + phi[i+1][j+1][k+1]*phi.dxc(i+1);
        real hy  = (hpc-hmc)/(dys(j)+dyn(j));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hyy = (hpc-2.0*hcc+hmc)/(phi.dyc(j)*phi.dyc(j));
        real hzz = (hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hyz = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dyc(j)*phi.dzc(k));
        //kappa[i][j][k] = copysign(1.0, nz)
        kappa[i][j][k] = -1.0
                       * (hyy + hzz + hyy*hz*hz + hzz*hy*hy + 2.0*hyz*hy*hz)
                       / pow(1.0 + hy*hy + hz*hz, 1.5);

      } else if (dirMax==3) {
        // calculate height
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
        hmm = phi[i-1][j-1][k-1]*phi.dzc(k-1)
            + phi[i-1][j-1][k  ]*phi.dzc(k)
            + phi[i-1][j-1][k+1]*phi.dzc(k+1);
        hcm = phi[i  ][j-1][k-1]*phi.dzc(k-1)
            + phi[i  ][j-1][k  ]*phi.dzc(k)
            + phi[i  ][j-1][k+1]*phi.dzc(k+1);
        hpm = phi[i+1][j-1][k-1]*phi.dzc(k-1)
            + phi[i+1][j-1][k  ]*phi.dzc(k)
            + phi[i+1][j-1][k+1]*phi.dzc(k+1);
        hmc = phi[i-1][j  ][k-1]*phi.dzc(k-1)
            + phi[i-1][j  ][k  ]*phi.dzc(k)
            + phi[i-1][j  ][k+1]*phi.dzc(k+1);
        hcc = phi[i  ][j  ][k-1]*phi.dzc(k-1)
            + phi[i  ][j  ][k  ]*phi.dzc(k)
            + phi[i  ][j  ][k+1]*phi.dzc(k+1);
        hpc = phi[i+1][j  ][k-1]*phi.dzc(k-1)
            + phi[i+1][j  ][k  ]*phi.dzc(k)
            + phi[i+1][j  ][k+1]*phi.dzc(k+1);
        hmp = phi[i-1][j+1][k-1]*phi.dzc(k-1)
            + phi[i-1][j+1][k  ]*phi.dzc(k)
            + phi[i-1][j+1][k+1]*phi.dzc(k+1);
        hcp = phi[i  ][j+1][k-1]*phi.dzc(k-1)
            + phi[i  ][j+1][k  ]*phi.dzc(k)
            + phi[i  ][j+1][k+1]*phi.dzc(k+1);
        hpp = phi[i+1][j+1][k-1]*phi.dzc(k-1)
            + phi[i+1][j+1][k  ]*phi.dzc(k)
            + phi[i+1][j+1][k+1]*phi.dzc(k+1);

        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hy  = (hcp-hcm)/(dys(j)+dyn(j));
        //real hxx =  2.0*h[-1][0]/(dxe(i)*(dxw(i)+dxe(i)))
        //          - 2.0*h[ 0][0]/(dxe(i)*dxw(i))
        //          + 2.0*h[-1][0]/(dxw(i)*(dxw(i)+dxe(i)));
        real hxx = (hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hyy = (hcp-2.0*hcc+hcm)/(phi.dyc(j)*phi.dyc(j));
        real hxy = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dxc(i)*phi.dyc(j));

        //kappa[i][j][k] = copysign(1.0, nz)
        kappa[i][j][k] = -1.0
                       * (hxx + hyy + hxx*hy*hy + hyy*hx*hx + 2.0*hxy*hx*hy)
                       / pow(1.0 + hx*hx + hy*hy, 1.5);
        //if(j==1&&(i==50||i==51)&&k==26){
        //if(j==1){
        //  std::cout<<i<<" "<<k<<" "<<kappa[i][j][k]<<" "<<copysign(1.0, nz)<<" "<<hx<<" "<<hxx<<"\n";;
        //}
      }
    }
  }
  kappa.exchange();
#endif

#if 0
  for_aijk(i,j,k) {
    stmp[i][j][k]=real(iflag[i][j][k]);
  }
  boil::plot->plot(phi,kappa,stmp, "phi-kappa-iflag", time->current_step());
  exit(0);
#endif

  return;
}
/*-----------------------------------------------------------------------------+
