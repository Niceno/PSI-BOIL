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
*                      No pre-coloring is considered.
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

/******************************************************************************/
void VOF::insert_bc_curv_HFnormal(const Scalar & scp, 
                                  const Comp cnormal,
                                  const Sign sig) {
/***************************************************************************//**
*  \brief Calculate curvature using height-function approach in 3D geometry.
*         Simple height extrapolations are used.
* 
*         Limitations: heights are constructed in the surface-normal direction
*                      (=> errors when CA ~ 90). Also, possible
*                      susceptibility to numerical errors creating spurious in-
*                      terfaces. Moreover, cells in second layer are not treated
*                      -> this is incorrect for high contact angles. 
*
*                      Wall assumed in negative z-dir.
*                      Detachment model is not considered.
*                      No pre-coloring is considered.
*
*     output: kappa
*******************************************************************************/

  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());

  arr2D heights, distances;
  heights.resize(hf_set.minorext);
  distances.resize(hf_set.minorext);
  for(auto & h : heights)
    h.resize(hf_set.minorext);
  for(auto & d : distances)
    d.resize(hf_set.minorext);
  real max_n; /* normal vector must be of good quality!!! */
  real mult;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
        bool test =  ((scp[i-1][j][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                   ||((scp[i-1][j][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                   ||((scp[i][j-1][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                   ||((scp[i][j+1][k]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                   ||((scp[i][j][k-1]-phisurf)*(scp[i][j][k]-phisurf)<=0.0)
                   ||((scp[i][j][k+1]-phisurf)*(scp[i][j][k]-phisurf)<=0.0);
        if(!test)
          continue;
        real dz0 = scp.dzc(k);
        for(auto & h : heights)
          for(auto & val : h)
            val = 0.0;
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
          for(int ii(-hf_set.nof); ii<=hf_set.nof; ++ii)
            for(int jj(-hf_set.nof); jj<=hf_set.nof; ++jj)
              heights[ii+hf_set.nof][jj+hf_set.nof] += dz * 
                (mult < 0 ? (1.-bounded_color(scp[i+ii][j+jj][k+kk]))
                          :     bounded_color(scp[i+ii][j+jj][k+kk]) );
        }

#ifdef USE_UPDATE_AT_WALLS
        for(auto & h : heights)
          for(auto & val : h)
            val -= dz0;
#endif
    
        real & hcc = heights[hf_set.nof][hf_set.nof];
        if(tol_wall*dz0<hcc&&hcc<=dz0) {
          tempflag[i][j][k] = 1;

          for(auto & h : heights)
            for(auto & val : h)
              if(val<tol_wall*dz0)
                val = 0.;

          for(int ii(-hf_set.nof); ii<=hf_set.nof; ++ii)
            for(int jj(-hf_set.nof); jj<=hf_set.nof; ++jj)
              distances[ii+hf_set.nof][jj+hf_set.nof] = 
                sqrt( (scp.xc(i+ii)-scp.xc(i))*(scp.xc(i+ii)-scp.xc(i))
                     +(scp.yc(j+jj)-scp.yc(j))*(scp.yc(j+jj)-scp.yc(j)) );
 
          kappa[i][j][k] = wall_curv_HFnormal_kernel(heights,
                                                     distances,
                                                     { scp.dxw(i),
                                                       scp.dxc(i),
                                                       scp.dxe(i) },
                                                     { scp.dys(j),
                                                       scp.dyc(j),
                                                       scp.dyn(j) },
                                                     ifull, jfull,
                                                     mult, max_n, cangle);

        } else {
          tempflag[i][j][k] = 0;
          kappa[i][j][k] = boil::unreal;
        }
      } /* next to off */
    } /* is on */
  } /* ijk */

  kappa.exchange();
  tempflag.exchange();

  return;
}
