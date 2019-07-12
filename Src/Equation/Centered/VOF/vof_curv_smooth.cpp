#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::curv_smooth() {
/***************************************************************************//**
*  \brief Calculate curvature using height function.
*     output: kappa
*******************************************************************************/

  smooth(phi,clr,4);

  gradphic(clr);
  
  set_iflag();

  /* div.norm */
  for_ijk(i,j,k) {
    if(abs(iflag[i][j][k])>nlayer-3){
    //if(abs(iflag[i][j][k])!=0){
      kappa[i][j][k]=boil::unreal;
    } else {
      real nw = nx[i-1][j][k];
      real ne = nx[i+1][j][k];
      real ns = ny[i][j-1][k];
      real nn = ny[i][j+1][k];
      real nb = nz[i][j][k-1];
      real nt = nz[i][j][k+1];

      kappa[i][j][k]=-((ne-nw)/(phi.dxw(i)+phi.dxe(i))
                      +(nn-ns)/(phi.dys(j)+phi.dyn(j))
                      +(nt-nb)/(phi.dzb(k)+phi.dzt(k)));
#if 0
      if(boil::cart.iam()==0&&i==9+2&&(j==1+2||j==2+2)&&k==16+2) {
        std::cout<<i<<" "<<j<<" "<<k<<" "<<kappa[i][j][k]<<"\n";
        std::cout<<ne<<" "<<nw<<" "
                 <<nt<<" "<<nb<<"\n";
      }
#endif
    }
  }

  kappa.exchange();
#ifdef DEBUG
  std::cout<<"curv::kappa "<<boil::cart.iam()<<"\n";
#endif

#if 0
  //extrapolate kappa
  for_aijk(i,j,k) {
    if(iflag[i][j][k]==0) {
      iflag[i][j][k]=1;
    } else {
      iflag[i][j][k]=0;
    }
  }
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
#if 0
  boil::plot->plot(kappa,nx,ny,nz, "kappa-nx-ny-nz", time->current_step());
  boil::plot->plot(sca,nx,ny,nz, "sca-nx-ny-nz", time->current_step());
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa", time->current_step());
  exit(0);
#endif
  /* curvature for symmetric plane (adjacent cells)
     exactly satisfy symmetricity */
  bnd_sym_kappa();
#ifdef DEBUG
  std::cout<<"curv::sym_kappa "<<boil::cart.iam()<<"\n";
#endif

  /* curvature for wall (adjacent cells and on boundary plane) 
     kappa in adjacent cell is necessary for interpolation, curv_interface() */
  bnd_wall_kappa();
#ifdef DEBUG
  std::cout<<"curv::bnd_wall_kappa "<<boil::cart.iam()<<"\n";
#endif

  /* boundary condition except wall (on boundary plane) */
  insert_bc_kappa(kappa);
  kappa.exchange_all();

#if 0
  boil::plot->plot(sca,nx,ny,nz, "sca-nx-ny-nz_beforeInter",
                     time->current_step());
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa_beforeInter",
                     time->current_step());
#endif

  /* interpolate curvature in interface cells */
  curv_interface();
#ifdef DEBUG
  std::cout<<"curv::curv_interface "<<boil::cart.iam()<<"\n";
#endif

#if 0
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa_afterInter",
                     time->current_step());
  boil::plot->plot(clr,iflag, "clr-iflag", time->current_step());
#endif

  /* extpolate curvature from interface cells to others */
  curv_interface_ext();

#if 0
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa_afterInterExt",
                     time->current_step());
  boil::plot->plot(clr,wflag, "clr-wflag", time->current_step());
  exit(0);
#endif
#endif

#if 0
  for_aijk(i,j,k) {
    stmp[i][j][k]=real(iflag[i][j][k]);
  }
  if(time->current_step() == 1) {
    boil::plot->plot(phi,kappa,stmp, "phi-kappa-iflag", time->current_step());
  }
  exit(0);
#endif
  return;
}

