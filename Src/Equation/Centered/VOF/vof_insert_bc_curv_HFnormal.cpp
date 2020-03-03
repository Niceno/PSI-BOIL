#include "vof.h"

#define USE_UPDATE_AT_WALLS

/******************************************************************************/
void VOF::insert_bc_curv_HFnormal(const Scalar & scp, 
                                  const Comp ctangential, const Comp cnormal,
                                  const Sign sig,
                                  const Range<int> ridx) {
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

  bool ranged(false);
  if(ridx.exists()) {
    ranged = true;
  }

  real hm,hc,hp;
  real max_n; /* normal vector must be of good quality!!! */
  real mult;
  for_i(i) {
    if(!ranged||ridx.contains(scp.domain()->global_I(i)-boil::BW+1)) {
     for_jk(j,k) {
       if(dom->ibody().on(i,j,k)) {
         if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
           bool test =  ((scp[i-1][j][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                      ||((scp[i-1][j][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                      ||((scp[i][j][k-1]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                      ||((scp[i][j][k+1]-phisurf)*(scp[i][j][k]-phisurf)<=0.0);
           if(!test)
             continue;
           real dz0 = scp.dzc(k);
           hm = hc = hp = 0.;
           max_n = -nz[i][j][k];
           if(max_n<0.0) {
             mult = -1;
           } else {
             mult =  1;
           }
#ifdef USE_UPDATE_AT_WALLS
           for(int kk(-1); kk<=hf_set.mof; ++kk) {
             real dz = scp.dzc(k+kk);
             if(kk==-1) dz = scp.dzc(k+kk+1);
#else
           for(int kk( 0); kk<=hf_set.mof; ++kk) {
             real dz = scp.dzc(k+kk);
#endif
             hm += (mult < 0 ? (1.-bounded_color(scp[i-1][j][k+kk]))
                             :     bounded_color(scp[i-1][j][k+kk]) ) * dz;
             hc += (mult < 0 ? (1.-bounded_color(scp[i  ][j][k+kk]))
                             :     bounded_color(scp[i  ][j][k+kk]) ) * dz;
             hp += (mult < 0 ? (1.-bounded_color(scp[i+1][j][k+kk]))
                             :     bounded_color(scp[i+1][j][k+kk]) ) * dz;
           }

#ifdef USE_UPDATE_AT_WALLS
           hm -= dz0;
           hc -= dz0;
           hp -= dz0;
#endif
        
           if(tol_wall*dz0<hc&&hc<=dz0) {
             tempflag[i][j][k] = 1;

             if(hm<tol_wall*dz0)
               hm = 0.;
             if(hp<tol_wall*dz0)
               hp = 0.;

             kappa[i][j][k] = wall_curv_HFnormal_kernel(scp.xc(i),
                                                        hm,hc,hp,
                                                        scp.dxw(i),
                                                        scp.dxc(i),scp.dxe(i),
                                                        mult, cangle);

             //boil::oout<<i<<" "<<j<<" "<<k<<" | "<<hm/dz0<<" "<<hc/dz0<<" "<<hp/dz0
             //          <<" | "<<mult<<" "<<max_n<<" "<<kappa[i][j][k]<<boil::endl;
           } else {
             tempflag[i][j][k] = 0;
             kappa[i][j][k] = boil::unreal;
           }
         }
       } /* is on */
     } /* jk */
    } /* in range */
  } /* i */

  kappa.exchange();
  tempflag.exchange();

  return;

}
