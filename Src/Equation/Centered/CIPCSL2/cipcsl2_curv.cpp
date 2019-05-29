#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void CIPCSL2::curv(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate curvature.
*     input : sca
*     output: kappa
*******************************************************************************/
  boil::timer.start("cipcsl2 curv");

  /* normal vector at cell center */
  gradphic(sca);

  /* div.norm */
  for_ijk(i,j,k) {
    if(abs(iflag[i][j][k])>nlayer-6){
      kappa[i][j][k]=0.0;
    } else {
      real nw = nx[i-1][j][k];
      real ne = nx[i+1][j][k];
      real ns = ny[i][j-1][k];
      real nn = ny[i][j+1][k];
      real nb = nz[i][j][k-1];
      real nt = nz[i][j][k+1];

      kappa[i][j][k]=-((ne-nw)/(sca.dxw(i)+sca.dxe(i))
                      +(nn-ns)/(sca.dys(j)+sca.dyn(j))
                      +(nt-nb)/(sca.dzb(k)+sca.dzt(k)));

    }
  }

#if 0
  boil::plot->plot(kappa,nx,ny,nz, "kappa-nx-ny-nz", time->current_step());
  boil::plot->plot(sca,nx,ny,nz, "sca-nx-ny-nz", time->current_step());
  boil::plot->plot(clr,sca,kappa, "clr-sca-kappa", time->current_step());
  exit(0);
#endif

  /* curvature for wall (adjacent cells and on boundary plane) 
     kappa in adjacent cell is necessary for interpolation, curv_interface() */
  bnd_wall_kappa();

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

  boil::timer.stop("cipcsl2 curv");
  return;
}

