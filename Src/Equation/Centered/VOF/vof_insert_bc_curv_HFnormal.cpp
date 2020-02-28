#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_curv_HFnormal(const Scalar & scp, 
                                  const Comp ctangential, const Comp cnormal,
                                  const Sign sig) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry.
* 
*         Reference: derived by me (Lubomir)
*
*         Limitations: heights are constructed in the surface-normal direction
*                      (=> errors when CA ~ 90). Also, possible
*                      susceptibility to numerical errors creating spurious in-
*                      terfaces. Moreover, cells in second layer are not treated
*                      -> this is incorrect for high contact angles. 
*
*                      Wall assumed in negative z-dir.
*                      Detachment model is not considered.
*
*
*     output: kappa
*******************************************************************************/

  assert(ctangential==Comp::i());
  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());

  real hm,hc,hp;
  real max_n; /* normal vector must be of good quality!!! */
  real mult;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
        real dz0 = scp.dzc(k);
        hm = hc = hp = 0.;
        max_n = -nz[i][j][k];
        if(max_n<0.0) {
          mult = -1;
        } else {
          mult =  1;
        }
        for(int kk( 0); kk<=hf_set.mof; ++kk) {
          real dz = scp.dzc(k+kk);
          if(kk==-1) dz = scp.dzc(k+kk+1);
          hm += (mult < 0 ? (1.-bounded_color(scp[i-1][j][k+kk]))
                          :     bounded_color(scp[i-1][j][k+kk]) ) * dz;
          hc += (mult < 0 ? (1.-bounded_color(scp[i  ][j][k+kk]))
                          :     bounded_color(scp[i  ][j][k+kk]) ) * dz;
          hp += (mult < 0 ? (1.-bounded_color(scp[i+1][j][k+kk]))
                          :     bounded_color(scp[i+1][j][k+kk]) ) * dz;
        }
        
        /* by subtracting dzc = ghost dzm, the threshold is zero */
        //hm -= dz0;
        //hc -= dz0;
        //hp -= dz0;

        if(0.0<hc&&hc<=dz0) {
          boil::oout<<i<<" "<<j<<" "<<k<<" | "<<hm<<" "<<hc<<" "<<hp<<" "<<mult_wall<<" "<<cangle<<boil::endl;

          tempflag[i][j][k] = 1;
          kappa[i][j][k] = wall_curv_HFnormal_kernel(scp.xc(i),
                                                     hm,hc,hp,
                                                     scp.dxw(i),
                                                     scp.dxc(i),scp.dxe(i),
                                                     mult_wall, cangle);
        } else {
          tempflag[i][j][k] = 0;
          kappa[i][j][k] = boil::unreal;
        }
      }
    } /* is on */
  } /* ijk */

  kappa.exchange();
  tempflag.exchange();

  return;

}
