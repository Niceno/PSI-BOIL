#include "vof.h"

#define USE_UPDATE_AT_WALLS

/******************************************************************************/
void VOF::insert_bc_curv_HFnormal(const Scalar & scp, 
                                  const Comp ctangential, const Comp cnormal,
                                  const Sign sig) {
                                  //const Range<int> ridx) {
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

#if 0
  bool ranged(false);
  if(ridx.exists()) {
    ranged = true;
  }
#endif

  real hm,hc,hp;
  real max_n; /* normal vector must be of good quality!!! */
  real mult;
  for_i(i) {
    if(true) {//!ranged||ridx.contains(scp.domain()->global_I(i)-boil::BW+1)) {
     for_jk(j,k) {
       if(dom->ibody().on(i,j,k)) {
         if(dom->ibody().off(i,j,k-1) || (k==sk() && bflag_struct.kminw)) {
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

             real cangval = cangle(i,j,k);
             kappa[i][j][k] = wall_curv_HFnormal_kernel(scp.xc(i),
                                                        hm,hc,hp,
                                                        scp.dxw(i),
                                                        scp.dxc(i),scp.dxe(i),
                                                        mult, cangval);

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
  real max_n,nnx,nny,dummy(0.); /* normal vector must be of good quality!!! */
  real mult;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      if(dom->ibody().off(i,j,k-1) || (k==sk() && bflag_struct.kminw)) {
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

        /* n pointing from the fluid which we are constructing heights of */
        max_n = -nz[i][j][k];
        if(max_n<0.0) {
          mult = -1;
          nnx  =  nx[i][j][k];
          nny  =  ny[i][j][k];
        } else {
          mult =  1;
          nnx  = -nx[i][j][k];
          nny  = -ny[i][j][k];
        }

        /* normalised projection on the x-y plane, dummy instead of z */
        normalize(nnx,nny,dummy);

        /* fill stencil */
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

        /* normalize stencil, so that zero appears at surface */
#ifdef USE_UPDATE_AT_WALLS
        for(auto & h : heights)
          for(auto & val : h)
            val -= dz0;
#endif
    
        /* we need to be above tolerance to consider contact with wall */
        real & hcc = heights[hf_set.nof][hf_set.nof];
        if(tol_wall*dz0<hcc&&hcc<=dz0) {
          tempflag[i][j][k] = 1;

          for(auto & h : heights)
            for(auto & val : h)
              if(val<tol_wall*dz0)
                val = 0.;

          /* distances for extrapolation of CA */
#if 0
          /* simple method, not taking into account x-y surface orientation,
             thus only using drop of height based on the contact-angle-based
             slope (CABS) */
          for(int ii(-hf_set.nof); ii<=hf_set.nof; ++ii)
            for(int jj(-hf_set.nof); jj<=hf_set.nof; ++jj)
              distances[ii+hf_set.nof][jj+hf_set.nof] = 
                sqrt( (scp.xc(i+ii)-scp.xc(i))*(scp.xc(i+ii)-scp.xc(i))
                     +(scp.yc(j+jj)-scp.yc(j))*(scp.yc(j+jj)-scp.yc(j)) );
#else
          /* a slope in any given direction should take into account the
             projection onto the boundary plane, calculating drop as 
                           CABS * < (nx,ny) , (Dx,Dy) >,
             where the latter term represents a scalar product of the normal
             vector projection onto the boundary plane and a distance vector
             in the chosen direction, n cdot D. Note that:
                           n cdot D = ||n|| * ||d|| * cos(alpha), 
              where ||n|| = 1 and alpha is the angle between the normal
              vector and the distance vector. Thus, n cdot D is numerically
              equal to the length of the projection of D in the n-dir */

          /* this is wrong. see curv_HF_kernel on how to implement */
          exit(0);
  
          /* here we calculate < (nx,ny) , (Dx,Dy) > and store it to be used
             by the kernel */
          for(int ii(-hf_set.nof); ii<=hf_set.nof; ++ii)
            for(int jj(-hf_set.nof); jj<=hf_set.nof; ++jj) {
              real distx = (scp.xc(i+ii)-scp.xc(i))*nnx;
              real disty = (scp.yc(j+jj)-scp.yc(j))*nny;

              /* we disallow extrapolation in the negative sense, this would
                 indicate corrupted normal vector */
              distances[ii+hf_set.nof][jj+hf_set.nof] = 
                std::max(0.,distx+disty);
            }
#endif
 
          real cangval = cangle(i,j,k);
          kappa[i][j][k] = wall_curv_HFnormal_kernel(heights,
                                                     distances,
                                                     { scp.dxw(i),
                                                       scp.dxc(i),
                                                       scp.dxe(i) },
                                                     { scp.dys(j),
                                                       scp.dyc(j),
                                                       scp.dyn(j) },
                                                     bflag_struct.ifull, bflag_struct.jfull,
                                                     mult, max_n, cangval);

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
