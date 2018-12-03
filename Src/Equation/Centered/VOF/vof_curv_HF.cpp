#include "vof.h"

/******************************************************************************/
void VOF::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  const int ist=-3;
  const int ied=+3;
  int kappa_flag =0;

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
  norm_cc(phi);
  // crude code: need boundary treatment for iflag //

  /* curvature calculation */
  kappa=0.0;

  for_aijk(i,j,k) {
    stmp[i][j][k] = phi[i][j][k];
  }

  for_ijk(i,j,k) {
    if (iflag[i][j][k]==1) {
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
    if(time->current_step() % 858 == 0) {
      std::cout<<"i: "<<i<<"j: "<<j<<"k: "<<k<<"dirMax: "<<dirMax<<"\n";     
    } 
#endif 
      
      if (dirMax==1) {
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;
#if 0
        hmm = phi[i-1][j-1][k-1]*phi.dxc(i-1)
            + phi[i  ][j-1][k-1]*phi.dxc(i)
            + phi[i+1][j-1][k-1]*phi.dxc(i+1);
#else
        hmm = 0.0, hcm = 0.0, hpm = 0.0;
        hmc = 0.0, hcc = 0.0, hpc = 0.0;
        hmp = 0.0, hcp = 0.0, hpp = 0.0;

        for (int ii=ist; ii<=-1; ii++) {  
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k-1]-phi[i+ii+1][j-1][k-1])>0){   //mm
            stmp[i+ii][j-1][k-1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k-1]-phi[i+ii+1][j][k-1])>0){   //cm
            stmp[i+ii][j][k-1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k-1]-phi[i+ii+1][j+1][k-1])>0){   //pm
            stmp[i+ii][j+1][k-1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k]-phi[i+ii+1][j-1][k])>0){   //mc
            stmp[i+ii][j-1][k]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k]-phi[i+ii+1][j][k])>0){   //cc
            stmp[i+ii][j][k]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k]-phi[i+ii+1][j+1][k])>0){   //pc
            stmp[i+ii][j+1][k]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k+1]-phi[i+ii+1][j-1][k+1])>0){   //mp
            stmp[i+ii][j-1][k+1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k+1]-phi[i+ii+1][j][k+1])>0){   //cp
            stmp[i+ii][j][k+1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k+1]-phi[i+ii+1][j+1][k+1])>0){   //pp
            stmp[i+ii][j+1][k+1]=0.5*(1+(-1)*copysign(1.0, nx[i][j][k]));
          }
        }

        for (int ii=1; ii<=ied; ii++) {
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k-1]-phi[i+ii-1][j-1][k-1])<0){   //mm
            stmp[i+ii][j-1][k-1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k-1]-phi[i+ii-1][j][k-1])<0){   //cm
            stmp[i+ii][j][k-1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k-1]-phi[i+ii-1][j+1][k-1])<0){   //pm
            stmp[i+ii][j+1][k-1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k]-phi[i+ii-1][j-1][k])<0){   //mc
            stmp[i+ii][j-1][k]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k]-phi[i+ii-1][j][k])<0){   //cc
            if (i==29 && j==20 && k==20 && ii==3){
              std::cout<<"stmp[32][20][20] = "<<stmp[32][20][20]<<"\n";
            }
            stmp[i+ii][j][k]=0.5*(1+copysign(1.0, nx[i][j][k]));
            if (i==29 && j==20 && k==20 && ii==3){
              std::cout<<"stmp[32][20][20] = "<<stmp[32][20][20]<<"\n";
            }  
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k]-phi[i+ii-1][j+1][k])<0){   //pc
            stmp[i+ii][j+1][k]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j-1][k+1]-phi[i+ii-1][j-1][k+1])<0){   //mp
            stmp[i+ii][j-1][k+1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j][k+1]-phi[i+ii-1][j][k+1])<0){   //cp
            stmp[i+ii][j][k+1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
          if (copysign(1.0, nx[i][j][k])*(phi[i+ii][j+1][k+1]-phi[i+ii-1][j+1][k+1])<0){   //pp
            stmp[i+ii][j+1][k+1]=0.5*(1+copysign(1.0, nx[i][j][k]));
          }
        }

        for (int ii=ist; ii<=ied; ii++) {
          hmm += stmp[i+ii][j-1][k-1]*phi.dxc(i+ii);
          hcm += stmp[i+ii][j][k-1]*phi.dxc(i+ii);
          hpm += stmp[i+ii][j+1][k-1]*phi.dxc(i+ii);
          hmc += stmp[i+ii][j-1][k]*phi.dxc(i+ii);
          hcc += stmp[i+ii][j][k]*phi.dxc(i+ii);
          hpc += stmp[i+ii][j+1][k]*phi.dxc(i+ii);
          hmp += stmp[i+ii][j-1][k+1]*phi.dxc(i+ii);
          hcp += stmp[i+ii][j][k+1]*phi.dxc(i+ii); 
          hpp += stmp[i+ii][j+1][k+1]*phi.dxc(i+ii);
        }
#endif
        if (ied*phi.dxc(i)<hcc && hcc<(ied+1)*phi.dxc(i)) {                  //Only for uniform mesh
#if 0  
          if (k=1) {
           std::cout<<"hcc "<<hcc<<" ied*phi.dxc(i) "<<ied*phi.dxc(i) <<\
           " (ied+1)*phi.dxc(i) "<<(ied+1)*phi.dxc(i)<<" i,j,k "<<i<<" "<<\
           j<<" "<<k<<"\n";
          }
#endif
        kappa_flag += 1;
        real hy  = (hpc-hmc)/(dys(j)+dyn(j));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));
        real hyy = (hpc-2.0*hcc+hmc)/(phi.dyc(j)*phi.dyc(j));
        real hzz = (hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hyz = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dyc(j)*phi.dzc(k));
        //kappa[i][j][k] = copysign(1.0, nz)
        kappa[i][j][k] = -1.0
                       * (hyy + hzz + hyy*hz*hz + hzz*hy*hy - 2.0*hyz*hy*hz)
                       / pow(1.0 + hy*hy + hz*hz, 1.5);
        }
#if 1
        if (i==29 && j==20 && k==20){
          std::cout<<"hcc "<<hcc<<"\n";
          std::cout<<"phi "<<phi[29][20][20]<<"\n";
          std::cout<<"kappa "<<kappa[29][20][20]<<"\n";   
          std::cout<<"nx "<<nx[29][20][20]<<" ny "<<ny[29][20][20]<<" nz "<<nz[29][20][20]<<"\n";
          std::cout<<"stmp[29][20][20] "<<stmp[29][20][20]<<"\n";
          std::cout<<"stmp[32][20][20] "<<stmp[32][20][20]<<"\n";
          std::cout<<"stmp[33][20][20] "<<stmp[33][20][20]<<"\n";
          std::cout<<"copysign(1.0, nx[i][j][k])*(stmp[i+ii][j][k]-stmp[i+ii-1][j][k])"<<     \
                     copysign(1.0, nx[29][20][20])*(stmp[32][20][20]-stmp[31][20][20])<<"\n";
        }
        if (i==30 && j==20 && k==20){
          std::cout<<"hcc "<<hcc<<"\n";
          std::cout<<"phi "<<phi[30][20][20]<<"\n";
          std::cout<<"kappa "<<kappa[30][20][20]<<"\n";
          std::cout<<"nx "<<nx[30][20][20]<<" ny "<<ny[30][20][20]<<" nz "<<nz[30][20][20]<<"\n";
        }

        if (i==12 && j==20 && k==20){
          std::cout<<"hcc "<<hcc<<"\n";
          std::cout<<"phi "<<phi[12][20][20]<<"\n";
          std::cout<<"kappa "<<kappa[12][20][20]<<"\n";
          std::cout<<"nx "<<nx[12][20][20]<<" ny "<<ny[12][20][20]<<" nz "<<nz[12][20][20]<<"\n";
        }
        if (i==11 && j==20 && k==20){
          std::cout<<"hcc "<<hcc<<"\n";
          std::cout<<"phi "<<phi[11][20][20]<<"\n";
          std::cout<<"kappa "<<kappa[11][20][20]<<"\n";
          std::cout<<"nx "<<nx[11][20][20]<<" ny "<<ny[11][20][20]<<" nz "<<nz[11][20][20]<<"\n";
        }


#endif

        for (int ii=ist; ii<=ied; ii++) {
          stmp[i+ii][j-1][k-1] = phi[i+ii][j-1][k-1];
          stmp[i+ii][j][k-1]   = phi[i+ii][j][k-1];
          stmp[i+ii][j+1][k-1] = phi[i+ii][j+1][k-1];
          stmp[i+ii][j-1][k]   = phi[i+ii][j-1][k];
          stmp[i+ii][j][k]     = phi[i+ii][j][k];
          stmp[i+ii][j+1][k]   = phi[i+ii][j+1][k];
          stmp[i+ii][j-1][k+1] = phi[i+ii][j-1][k+1];
          stmp[i+ii][j][k+1]   = phi[i+ii][j][k+1];
          stmp[i+ii][j+1][k+1] = phi[i+ii][j+1][k+1]; 
        } 

      } else if (dirMax==3) {
        // calculate height
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;

        hmm = 0.0, hcm = 0.0, hpm = 0.0;
        hmc = 0.0, hcc = 0.0, hpc = 0.0;
        hmp = 0.0, hcp = 0.0, hpp = 0.0;

        for (int ii=ist; ii<=-1; ii++) {
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j-1][k+ii]-phi[i-1][j-1][k+ii+1])>0){   //mm
            stmp[i-1][j-1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j-1][k+ii]-phi[i][j-1][k+ii+1])>0){   //cm
            stmp[i][j-1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j-1][k+ii]-phi[i+1][j-1][k+ii+1])>0){   //pm
            stmp[i+1][j-1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j][k+ii]-phi[i-1][j][k+ii+1])>0){   //mc
            stmp[i-1][j][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j][k+ii]-phi[i][j][k+ii+1])>0){   //cc
            stmp[i][j][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j][k+ii]-phi[i+1][j][k+ii+1])>0){   //pc
            stmp[i+1][j][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j+1][k+ii]-phi[i-1][j+1][k+ii+1])>0){   //mp
            stmp[i-1][j+1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j+1][k+ii]-phi[i][j+1][k+ii+1])>0){   //cp
            stmp[i][j+1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j+1][k+ii]-phi[i+1][j+1][k+ii+1])>0){   //pp
            stmp[i+1][j+1][k+ii]=0.5*(1+(-1)*copysign(1.0, nz[i][j][k]));
          }
        }

        for (int ii=1; ii<=ied; ii++) { 
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j-1][k+ii]-phi[i-1][j-1][k+ii-1])<0){   //mm
            stmp[i-1][j-1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j-1][k+ii]-phi[i][j-1][k+ii-1])<0){   //cm
            stmp[i][j-1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j-1][k+ii]-phi[i+1][j-1][k+ii-1])<0){   //pm
            stmp[i+1][j-1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j][k+ii]-phi[i-1][j][k+ii-1])<0){   //mc
            stmp[i-1][j][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j][k+ii]-phi[i][j][k+ii-1])<0){   //cc
            stmp[i][j][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j][k+ii]-phi[i+1][j][k+ii-1])<0){   //pc
            stmp[i+1][j][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i-1][j+1][k+ii]-phi[i-1][j+1][k+ii-1])<0){   //mp
            stmp[i-1][j+1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i][j+1][k+ii]-phi[i][j+1][k+ii-1])<0){   //cp
            stmp[i][j+1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
          if (copysign(1.0, nz[i][j][k])*(phi[i+1][j+1][k+ii]-phi[i+1][j+1][k+ii-1])<0){   //pp
            stmp[i+1][j+1][k+ii]=0.5*(1+copysign(1.0, nz[i][j][k]));
          }
        }  



        for (int ii=ist; ii<=ied; ii++) {    
          hmm += stmp[i-1][j-1][k+ii]*phi.dxc(i+ii);
          hcm += stmp[i][j-1][k+ii]*phi.dxc(i+ii);
          hpm += stmp[i+1][j-1][k+ii]*phi.dxc(i+ii);
          hmc += stmp[i-1][j][k+ii]*phi.dxc(i+ii);
          hcc += stmp[i][j][k+ii]*phi.dxc(i+ii);
          hpc += stmp[i+1][j][k+ii]*phi.dxc(i+ii);
          hmp += stmp[i-1][j+1][k+ii]*phi.dxc(i+ii); 
          hcp += stmp[i][j+1][k+ii]*phi.dxc(i+ii);
          hpp += stmp[i+1][j+1][k+ii]*phi.dxc(i+ii);
        }
        
        if (ied*phi.dxc(i)<hcc && hcc<(ied+1)*phi.dxc(i)) {
        kappa_flag += 1;
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
                       * (hxx + hyy + hxx*hy*hy + hyy*hx*hx - 2.0*hxy*hx*hy)
                       / pow(1.0 + hx*hx + hy*hy, 1.5);
        //if(j==1&&(i==50||i==51)&&k==26){
        //if(j==1){
        //  std::cout<<i<<" "<<k<<" "<<kappa[i][j][k]<<" "<<copysign(1.0, nz)<<" "<<hx<<" "<<hxx<<"\n";;
        //}
        }
        for (int ii=ist; ii<=ied; ii++) {
          stmp[i-1][j-1][k+ii] = phi[i-1][j-1][k+ii];
          stmp[i][j-1][k+ii]   = phi[i][j-1][k+ii];
          stmp[i+1][j-1][k+ii] = phi[i+1][j-1][k+ii];
          stmp[i-1][j][k+ii]   = phi[i-1][j][k+ii];
          stmp[i][j][k+ii]     = phi[i][j][k+ii];
          stmp[i+1][j][k+ii]   = phi[i+1][j][k+ii];
          stmp[i-1][j+1][k+ii] = phi[i-1][j+1][k+ii];
          stmp[i][j+1][k+ii]   = phi[i][j+1][k+ii];
          stmp[i+1][j+1][k+ii] = phi[i+1][j+1][k+ii];
        }
      } else if(dirMax==2){
        // calculate height
        real hmm, hcm, hpm, hmc, hcc, hpc, hmp, hcp, hpp;

        hmm = 0.0, hcm = 0.0, hpm = 0.0;
        hmc = 0.0, hcc = 0.0, hpc = 0.0;
        hmp = 0.0, hcp = 0.0, hpp = 0.0;

        for (int ii=ist; ii<=-1; ii++) { 
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k-1]-phi[i-1][j+ii+1][k-1])>0){  //mm
            stmp[i-1][j+ii][k-1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k-1]-phi[i][j+ii+1][k-1])>0){  //cm
            stmp[i][j+ii][k-1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k-1]-phi[i+1][j+ii+1][k-1])>0){  //pm
            stmp[i+1][j+ii][k-1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k]-phi[i-1][j+ii+1][k])>0){  //mc
            stmp[i-1][j+ii][k]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k]-phi[i][j+ii+1][k])>0){  //cc
            stmp[i][j+ii][k]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k]-phi[i+1][j+ii+1][k])>0){  //pc
            stmp[i+1][j+ii][k]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k+1]-phi[i-1][j+ii+1][k+1])>0){  //mp
            stmp[i-1][j+ii][k+1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k+1]-phi[i][j+ii+1][k+1])>0){  //cp
            stmp[i][j+ii][k+1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k+1]-phi[i+1][j+ii+1][k+1])>0){  //pp
            stmp[i+1][j+ii][k+1]=0.5*(1+(-1)*copysign(1.0, ny[i][j][k]));
          }
        }
        
        for (int ii=1; ii<=ied; ii++) {
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k-1]-phi[i-1][j+ii-1][k-1])<0){   //mm
            stmp[i-1][j+ii][k-1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k-1]-phi[i][j+ii-1][k-1])<0){   //cm
            stmp[i][j+ii][k-1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k-1]-phi[i+1][j+ii-1][k-1])<0){   //pm
            stmp[i+1][j+ii][k-1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k]-phi[i-1][j+ii-1][k])<0){   //mc
            stmp[i-1][j+ii][k]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k]-phi[i][j+ii-1][k])<0){   //cc
            stmp[i][j+ii][k]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k]-phi[i+1][j+ii-1][k])<0){   //pc
            stmp[i+1][j+ii][k]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i-1][j+ii][k+1]-phi[i-1][j+ii-1][k+1])<0){   //mp
            stmp[i-1][j+ii][k+1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i][j+ii][k+1]-phi[i][j+ii-1][k+1])<0){   //cp
            stmp[i][j+ii][k+1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
          if (copysign(1.0, ny[i][j][k])*(phi[i+1][j+ii][k+1]-phi[i+1][j+ii-1][k+1])<0){   //pp
            stmp[i+1][j+ii][k+1]=0.5*(1+copysign(1.0, ny[i][j][k]));
          }
        }

        for (int ii=ist; ii<=ied; ii++) {
          hmm += stmp[i-1][j+ii][k-1]*phi.dxc(i+ii);
          hcm += stmp[i][j+ii][k-1]*phi.dxc(i+ii);
          hpm += stmp[i+1][j+ii][k-1]*phi.dxc(i+ii);
          hmc += stmp[i-1][j+ii][k]*phi.dxc(i+ii);
          hcc += stmp[i][j+ii][k]*phi.dxc(i+ii);
          hpc += stmp[i+1][j+ii][k]*phi.dxc(i+ii);
          hmp += stmp[i-1][j+ii][k+1]*phi.dxc(i+ii);
          hcp += stmp[i][j+ii][k+1]*phi.dxc(i+ii);
          hpp += stmp[i+1][j+ii][k+1]*phi.dxc(i+ii);
        }

        if (ied*phi.dxc(i)<hcc && hcc<(ied+1)*phi.dxc(i)) {
        
        kappa_flag += 1;
        real hx  = (hpc-hmc)/(dxw(i)+dxe(i));
        real hz  = (hcp-hcm)/(dzb(k)+dzt(k));

        real hxx = (hpc-2.0*hcc+hmc)/(phi.dxc(i)*phi.dxc(i));
        real hzz = (hcp-2.0*hcc+hcm)/(phi.dzc(k)*phi.dzc(k));
        real hxz = (hpp-hpm-hmp+hmm)
                 / (4.0*phi.dxc(i)*phi.dzc(k));
        kappa[i][j][k] = -1.0
                       * (hxx + hzz + hxx*hz*hz + hzz*hx*hx - 2.0*hxz*hx*hz)
                       / pow(1.0 + hx*hx + hz*hz, 1.5);                      
        }
        for (int ii=ist; ii<=ied; ii++) {
          stmp[i-1][j+ii][k-1] = phi[i-1][j+ii][k-1];
          stmp[i][j+ii][k-1]   = phi[i][j+ii][k-1];
          stmp[i+1][j+ii][k-1] = phi[i+1][j+ii][k-1];
          stmp[i-1][j+ii][k]   = phi[i-1][j+ii][k];
          stmp[i][j+ii][k]     = phi[i][j+ii][k];
          stmp[i+1][j+ii][k]   = phi[i+1][j+ii][k];
          stmp[i-1][j+ii][k+1] = phi[i-1][j+ii][k+1];
          stmp[i][j+ii][k+1]   = phi[i][j+ii][k+1];
          stmp[i+1][j+ii][k+1] = phi[i+1][j+ii][k+1];
        }
      }
    }
  }
  kappa.exchange();

#if 1
  for_aijk(i,j,k) {
    stmp[i][j][k]=real(iflag[i][j][k]);
  }
  if(time->current_step() == 1) {
    boil::plot->plot(phi,kappa,stmp, "phi-kappa-iflag", time->current_step());
  }
//  exit(0);
#endif

  //calculate maximum kappa and averaged kappa for one sphere or circle
#if 0
  real kappa_max = 0.0;
  real kappa_exact = 1.0;
  real L1    = 0.0;
  real L2    = 0.0;
  real L_inf = 0.0;
  real kappa_average;
  for_ijk(i,j,k) {
    if (kappa[i][j][k]>kappa_max) {
      kappa_max = kappa[i][j][k];  
    }
    if (kappa[i][j][k]>(kappa_exact/2.0)) {       
      L1 += fabs(kappa[i][j][k]-kappa_exact);                 
      L2 += pow((kappa[i][j][k] - kappa_exact),2.0);
    }
    if(kappa[i][j][k]>(kappa_exact/2.0)) {
      if (fabs(kappa[i][j][k]-kappa_exact)>L_inf) {
        L_inf = fabs(kappa[i][j][k]-kappa_exact);  
      }
    }
    kappa_average += kappa[i][j][k];
  }
  kappa_average = kappa_average/kappa_flag;
  L1            = L1/kappa_flag;
  L2            = sqrt(L2/kappa_flag);
  std::cout<<"kappa_max "<<kappa_max<<"\n";
  std::cout<<"kappa_average "<<kappa_average<<" kappa_flag "<<kappa_flag<<"\n";
  std::cout<<"L1 "<<L1<<"\n";
  std::cout<<"L2 "<<L2<<"\n";
  std::cout<<"L_inf "<<L_inf<<"\n";
#endif

  //calculate maximum kappa and averaged kappa for donut
#if 1
  kappa_flag = 0;
  real kappa_exact = 1.0;
  real L1         = 0.0;
  real L_inf = 0.0;
  for_ijk(i,j,k){
    if(kappa[i][j][k]>(kappa_exact/2.0)) {
      kappa_flag += 1;
    }
    if (kappa[i][j][k]>(kappa_exact/2.0)) {
      L1 += fabs(kappa[i][j][k]-kappa_exact);
    }
    if(kappa[i][j][k]>(kappa_exact/2.0)) {
      if (fabs(kappa[i][j][k]-kappa_exact)>L_inf) {
        L_inf = fabs(kappa[i][j][k]-kappa_exact);
      }
    }
  }
  L1            = L1/kappa_flag;
  std::cout<<"Calculate maximum kappa and averaged kappa for donut"<<"\n";
  std::cout<<"L1 "<<L1<<"\n";  
  std::cout<<"L_inf "<<L_inf<<"\n";
  std::cout<<"kappa_flag "<<kappa_flag<<"\n"; 
#endif  

#if 1
  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();
  if(time->current_step() == 1) {
    boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  }
#endif
  return;
}

