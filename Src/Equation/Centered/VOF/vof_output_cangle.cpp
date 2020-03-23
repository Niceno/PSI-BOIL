#include "vof.h"

/******************************************************************************/
void VOF::output_cangle_2d(const Comp ctangential, const Comp cnormal,
                           const Sign sig) {

  real h0(0.0), h1(0.0), h2(0.0);
  real dzzt0(0.0), dzzc0(0.0), dzzt1(0.0), dzzc1(0.0); 

  Range<int> ridx;

  output_cangle_2d(color(),ctangential,cnormal,
                   sig,ridx,h0,h1,h2,
                   dzzt0,dzzc0,dzzt1,dzzc1);

  return;

}

/******************************************************************************/
void VOF::output_cangle_2d(const Scalar & scp, 
                           const Comp ctangential, const Comp cnormal,
                           const Sign sig, const Range<int> ridx,
                           real & h0, real & h1, real & h2,
                           real & dzzt0, real & dzzc0,
                           real & dzzt1, real & dzzc1) {
/***************************************************************************//**
*  \brief Output contact angle for a 2D system. 
*         Also used by insertbc_curvHF_parallel.
*
*         Limitations: only single bubble/droplet with no overlap of near-wall
*                      cells is assumed (= no significant inward bending of 
*                      interface), wall assumed in negative z-dir.
*
*******************************************************************************/

  assert(ctangential==Comp::i());
  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());

  bool ranged(false);
  if(ridx.exists()) {
    ranged = true;
  }

  h0 = h1 = h2 = dzzt0 = dzzc0 = dzzt1 = dzzc1 = 0.0;
  for_i(i) {
    if(!ranged||ridx.contains(scp.domain()->global_I(i)-boil::BW+1)) {
      for_jk(j,k) {
        if(dom->ibody().on(i,j,k)) {
          if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
            h0 += (mult_wall < 0 ? (1.-scp[i][j][k  ]) : scp[i][j][k  ]) * scp.dxc(i);
            h1 += (mult_wall < 0 ? (1.-scp[i][j][k+1]) : scp[i][j][k+1]) * scp.dxc(i);
            h2 += (mult_wall < 0 ? (1.-scp[i][j][k+2]) : scp[i][j][k+2]) * scp.dxc(i);
            /* this should be the same for all cells! */
            dzzt0 = scp.dzt(k);
            dzzc0 = scp.dzc(k);
            dzzt1 = scp.dzt(k+1);
            dzzc1 = scp.dzc(k+1);
          }
        } /* is on */
      } /* jk */
    } /* in range */
  } /* i */
  boil::cart.sum_real(&h0);
  boil::cart.sum_real(&h1);
  boil::cart.sum_real(&h2);
  boil::cart.max_real(&dzzt0);
  boil::cart.max_real(&dzzc0);
  boil::cart.max_real(&dzzt1);
  boil::cart.max_real(&dzzc1);

  boil::oout<<"VOF::bccurv-h0: "<<time->current_time()<<" "<<h0<<" "<<h1<<" "<<180./boil::pi*atan(fabs(h1-h0)/dzzt0)<<boil::endl;

  return;
}
