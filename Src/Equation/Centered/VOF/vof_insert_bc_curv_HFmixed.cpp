#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_curv_HFmixed(const Scalar & scp, 
                                 const Comp ctangential, const Comp cnormal,
                                 const Sign sig) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry. The detachment criterion is set.
* 
*         Reference: derived by me (Lubomir)
*
*         Limitations: only single bubble/droplet with no overlap of near-wall
*                      cells is assumed (= no significant inward bending of 
*                      interface), wall assumed in negative z-dir.
*
*         This is resulting from the "infinity x 2" stencil used.
*
*     output: kappa
*******************************************************************************/

  assert(ctangential==Comp::i());
  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());

  //int errcnt(0);
  real h0(0.0), h1(0.0), h2(0.0);
  real dzzt0(0.0), dzzc0(0.0), dzzt1(0.0), dzzc1(0.0);
  for_ijk(i,j,k) {
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
  } /* ijk */
  boil::cart.sum_real(&h0);
  boil::cart.sum_real(&h1);
  boil::cart.sum_real(&h2);
  boil::cart.max_real(&dzzt0);
  boil::cart.max_real(&dzzc0);
  boil::cart.max_real(&dzzt1);
  boil::cart.max_real(&dzzc1);

  boil::oout<<"VOF::bccurv-h0: "<<time->current_time()<<" "<<h0<<" "<<h1<<" "<<180./boil::pi*atan(fabs(h1-h0)/dzzt0)<<boil::endl;

  if(!detachment_model.initialized()||!detachment_model.detached()) {

    real kappa_wall = wall_curv_HFmixed_kernel(h0,h1,dzzc0,dzzt0,
                                               mult_wall,
                                               cangle);
    
    real kappa_above = wall_curv_HFmixed_kernel(h0,h1,h2,
                                                dzzt0,dzzc1,dzzt1,
                                                mult_wall);
   
    for_ijk(i,j,k) {
      if(dom->ibody().on(i,j,k)) {
        if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
          if( (scp.xn(i  )<h0) && (h0<scp.xn(i+1)) ) {
            kappa[i][j][k] = kappa_wall;
            tempflag[i][j][k] = 1;
#if 0
          } else {
            kappa[i][j][k] = kappa_wall;
            tempflag[i][j][k] = 2;
#endif
          }
          if( (scp.xn(i  )<h1) && (h1<scp.xn(i+1)) ) {
            kappa[i][j][k+1] = kappa_above;
            tempflag[i][j][k+1] = 1;
#if 0
          } else {
            kappa[i][j][k+1] = kappa_above;
            tempflag[i][j][k+1] = 2;
#endif
          }
        } /* next to wall */
      } /* is on */
    } /* ijk */

  } /* ?detached */

  /* detachment status is updated for the next time step */
  if(detachment_model.initialized())
    detachment_model.test_detachment(h0/dzzc0);

  kappa.exchange();
  tempflag.exchange();

  return;

}
