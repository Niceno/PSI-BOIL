#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_curv_HFparallel(const Scalar & scp, 
                                    const Comp ctangential, const Comp cnormal,
                                    const Sign sig,
                                    const Range<int> ridx) {
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

  bool ranged(false);
  if(ridx.exists()) {
    ranged = true;
  }

  real h0(0.0), h1(0.0), h2(0.0);
  real dzzt0(0.0), dzzc0(0.0), dzzt1(0.0), dzzc1(0.0);

  output_cangle_2d(scp,ctangential,cnormal,
                   sig,ridx,h0,h1,h2,
                   dzzt0,dzzc0,dzzt1,dzzc1);

  if(!detachment_model.initialized()||!detachment_model.detached()) {

    real kappa_wall = wall_curv_HFparallel_kernel(h0,h1,dzzc0,dzzt0,
                                                  mult_wall,
                                                  cangle);
    
    real kappa_above = wall_curv_HFparallel_kernel(h0,h1,h2,
                                                   dzzt0,dzzc1,dzzt1,
                                                   mult_wall);
   
    for_i(i) {
      if(!ranged||ridx.contains(scp.domain()->global_I(i)-boil::BW+1)) {
        for_jk(j,k) {
          if(dom->ibody().on(i,j,k)) {
            if(dom->ibody().off(i,j,k-1) || (k==sk() && bflag_struct.kminw)) {
              if( (scp.xn(i  )<h0) && (h0<scp.xn(i+1)) ) {
                kappa[i][j][k] = kappa_wall;
                tempflag[i][j][k] = 1;
#if 1
              } else {
                kappa[i][j][k] = kappa_wall;
                tempflag[i][j][k] = 2;
#else
              } else if( (scp.xn(i+1)<h0) && (h0<scp.xn(i+2)) ) {
                kappa[i][j][k] = kappa_wall;
                tempflag[i][j][k] = 2;
              } else if( (scp.xn(i-1)<h0) && (h0<scp.xn(i  )) ) {
                kappa[i][j][k] = kappa_wall;
                tempflag[i][j][k] = 2;
#endif
              }
              if( (scp.xn(i  )<h1) && (h1<scp.xn(i+1)) ) {
                kappa[i][j][k+1] = kappa_above;
                tempflag[i][j][k+1] = 1;
#if 1
              } else {
                kappa[i][j][k+1] = kappa_above;
                tempflag[i][j][k+1] = 2;
#endif
              }
            } /* next to wall */
          } /* is on */
        } /* jk */
      } /* in range */
    } /* i */

  } /* ?detached */

  /* detachment status is updated for the next time step */
  if(detachment_model.initialized())
    detachment_model.test_detachment(h0/dzzc0);

  kappa.exchange();
  tempflag.exchange();

  return;

}
