#include "vofaxisym.h"
#define LOPEZ

/******************************************************************************/
void VOFaxisym::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function in 2D axisymmetric.
*     7 x 3 stencil
*     J.Lopez et al., Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*
*     modification: selection of cells where curvature is calculated is 
*                   modified from Lopez et al.
*
*     modification: adapted for a 2D axisymmetric grid
*
*     output: kappa
*******************************************************************************/

  /*-----------------------------------------------------------+
  |  Step 1: calculate normal vector                           |
  |  Step 2: define iflag=1                                    |
  |  Step 3: calculate curvature if (iflag==1 & 3<height<=4)   |
  |  Step 4: calculate curvature near wall (iflag = 10)        |
  |  Step 5: extrapolate curvature to five layers              |
  +-----------------------------------------------------------*/

  /* prepare 3x7 stencils: one for heights, one for dx */
  arr stencil, gridstencil;
  stencil.resize(3);
  gridstencil.resize(3);
  for(auto & s : stencil) s.resize(7);
  for(auto & g : gridstencil) g.resize(7);

  /* Normal vector is used for (i) mMax and (ii) Eq. (2) in J.Lopez et al. */
  if(norm_method_curvature != norm_method_advance) {
    norm(clr,norm_method_curvature,false); /* alpha is not extracted */
  }

  /* nx, ny, nz themselves are changed */
  /* mx, my, mz remain unchanged */
  true_norm_vect(nx,ny,nz,nx,ny,nz);

  /* define iflag */
  iflag=0;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      /* exclude near-wall cells: treated specially */
      if(   (i==si() && iminw) || (i==ei() && imaxw)
         || (k==sk() && kminw) || (k==ek() && kmaxw)
         || dom->ibody().off(i-1,j,k) || dom->ibody().off(i+1,j,k)
         || dom->ibody().off(i,j,k-1) || dom->ibody().off(i,j,k+1)
        ) {
        iflag[i][j][k]=10;
      } else {
        if((clr[i-1][j][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ iflag[i][j][k]=1;}
        if((clr[i+1][j][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ iflag[i][j][k]=1;}
        if((clr[i][j][k-1]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ iflag[i][j][k]=1;}
        if((clr[i][j][k+1]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ iflag[i][j][k]=1;}
      }
    }
  }
  iflag.bnd_update_symmetry(); /* copy on symmetry plane */
  iflag.exchange();

  /*------------------------+
  |  curvature calculation  |
  +------------------------*/
  kappa=boil::unreal;

  for_ijk(i,j,k) {

    /* calculate curvature only when iflag=1 */
    if (iflag[i][j][k]==1) {

      /* select direction of stencil-7 */
      Comp mMax;
      /* n points to the liquid */
      real nnx = -nx[i][j][k];
      real nnz = -nz[i][j][k];
      real abs_nx = fabs(nnx);
      real abs_nz = fabs(nnz);
      real max_abs_n, max_n;

      if (abs_nx<abs_nz) {
        mMax = Comp::k();
        max_abs_n = abs_nz;
        max_n = nnz;
      } else {
        mMax = Comp::i();
        max_abs_n = abs_nx;
        max_n = nnx;
      }

      if (mMax==Comp::i()) {
        if(!ifull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int imin=-3;
        int imax=3;  /* normal stencil size 7 (=3+1+3) */

        /* limit stencil size for cut-stencil */
        if(iminc) imin=std::max(-3,si()-i);
        if(imaxc) imax=std::min( 3,ei()-i);

        if(dom->ibody().off(i-2,j,k)) { imin=-1; }
        else if(dom->ibody().off(i-3,j,k)) { imin=-2; }

        if(dom->ibody().off(i+2,j,k)) { imax=1; }
        else if(dom->ibody().off(i+3,j,k)) { imax=2; }

        /* fill stencil */
        for(int kk(-1); kk<=1; ++kk) {
          /* stencil is reset with negative values representing undefined vals */
          for(auto & s : stencil[kk+1]) s = -1.0;
          for(auto & g : gridstencil[kk+1]) g = 0.0;
          for(int ii(imin); ii<=imax; ++ii) {
            gridstencil[kk+1][ii+3] = clr.dxc(i+ii);

            stencil[kk+1][ii+3] = std::min(1.0,std::max(0.0,clr[i+ii][j][k+kk]));
          }
        }

        /* calculate curvature */
        curv_HF_kernel_axisymmetric(stencil,gridstencil,
                                    imin,imax,
                                    clr.dzb(k),clr.dzc(k),clr.dzt(k),
                                    max_n,kfull,
                                    kappa[i][j][k],iflag[i][j][k]);

      } else {
        if(!kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int kmin=-3;
        int kmax=3;  /* normal stencil size 7 (=3+1+3) */

        /* limit stencil size for cut-stencil */
        if(kminc) kmin=std::max(-3,sk()-k);
        if(kmaxc) kmax=std::min( 3,ek()-k);

        if(dom->ibody().off(i,j,k-2)) { kmin=-1; }
        else if(dom->ibody().off(i,j,k-3)) { kmin=-2; }

        if(dom->ibody().off(i,j,k+2)) { kmax=1; }
        else if(dom->ibody().off(i,j,k+3)) { kmax=2; }

        /* fill stencil */
        for(int ii(-1); ii<=1; ++ii) {
          /* stencil is reset with negative values representing undefined vals */
          for(auto & s : stencil[ii+1]) s = -1.0;
          for(auto & g : gridstencil[ii+1]) g = 0.0;
          for(int kk(kmin); kk<=kmax; ++kk) {
            gridstencil[ii+1][kk+3] = clr.dzc(k+kk);

            stencil[ii+1][kk+3] = std::min(1.0,std::max(0.0,clr[i+ii][j][k+kk]));
          }
        }

        /* calculate curvature */
        curv_HF_kernel_axisymmetric(stencil,gridstencil,
                                    kmin,kmax,
                                    clr.dxw(i),clr.dxc(i),clr.dxe(i),
                                    max_n,ifull,
                                    kappa[i][j][k],iflag[i][j][k]);
      }  

    } /* iflag = 1 */
  } /* for ijk */

  kappa.bnd_update_symmetry(); /* copy on symmetry plane */
  iflag.bnd_update_symmetry(); /* copy on symmetry plane */
  kappa.exchange();
  iflag.exchange();

  /* near-wall calculations */
  if(!use_HF_wall) {
    boil::oout<<"VOFaxisym::curv_HF: Obsolete, underdeveloped code. Exiting."
              <<boil::endl;
    exit(0);
    //bdcurv();
  } else {

  }

  /* extrapolate kappa */
  stmp  = kappa;
  jflag = iflag;
  
  for(int iloop=1; iloop<5; iloop++) { /* 2019.07.09 */
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(iflag[i][j][k]==0) {
        /* at this point, the near-wall iflag must be properly 0 or 1 */
        int inb = std::min(1,iflag[i-1][j][k]) + std::min(1,iflag[i+1][j][k])
                + std::min(1,iflag[i][j][k-1]) + std::min(1,iflag[i][j][k+1]);
        if(inb >= 1) {
          stmp[i][j][k] = (real(std::min(1,iflag[i-1][j][k])) * kappa[i-1][j][k]
                        +  real(std::min(1,iflag[i+1][j][k])) * kappa[i+1][j][k]
                        +  real(std::min(1,iflag[i][j][k-1])) * kappa[i][j][k-1]
                        +  real(std::min(1,iflag[i][j][k+1])) * kappa[i][j][k+1])
                        /real(inb);
          jflag[i][j][k] = 2;  /* iflag=2 for extrapolated */
        }
      }
    }
    stmp.bnd_update_symmetry(); /* copy on symmetry plane */
    jflag.bnd_update_symmetry(); /* copy on symmetry plane */
    stmp.exchange();
    jflag.exchange();
    kappa = stmp;
    iflag = jflag;
  }
  
  return;
}

/*-------------------------------+
|  curvature calculation kernel  |
+-------------------------------*/
/* inputs:
   - stencil of color values
   - stencil of grid spacings in the dominant direction
   - starting and ending indices in the dominant direction
   - grid spacings in the minor direction
   - normal vector component in the dominant direction, pointing to gas
   - bool indicator if minor direction is a non-pseudo direction
   
   outputs:
   - kappa value (if computed)
   - iflag value (if kappa not computed)
*/
void VOFaxisym::curv_HF_kernel_axisymmetric(
                               arr & stencil, const arr & gridstencil,
                               const int imin, const int imax,
                               const real dm, const real dc, const real dp,
                               const real max_n, const bool truedir,
                               real & kap, int & flag
                                           ) {

  /* color is inverted if normal vector points in negative dir */
  real mult(1.0);
  if(max_n<0.0) {
    mult = -1.0;
    for(auto & s : stencil) {
      for(int ii(imin); ii<=imax; ++ii) {
        s[ii+3] = -s[ii+3];
      }
    }
  }

  /* Calculate Eq. (2) in J.Lopez et al. */
  /* Warning! if max_n is low accuracy (wrong sign), this algorithm
     doesn't work !!!  (Yohei) */
#ifdef LOPEZ
  for(auto & s : stencil) {
    for(int ii(-1); ii>=imin; ++ii) {
      if( (s[ii+3]-s[ii+4]) > 0.0) {
        if(  (s[ii+3]-phisurf)*(s[ii+4]-phisurf) > 0.0
           &&(s[ii+3]-phisurf)*(s[3   ]-phisurf) > 0.0) {
          s[ii+3] = s[ii+4]; /* Yohei */
        } else {
          s[ii+3] = 1.0; /* original */
        }
      }
    }

    for(int ii( 1); ii<=imax; ++ii) {
      if( (s[ii+3]-s[ii+2]) > 0.0) {
        if(  (s[ii+3]-phisurf)*(s[ii+2]-phisurf) > 0.0
           &&(s[ii+3]-phisurf)*(s[3   ]-phisurf) > 0.0) {
          s[ii+3] = s[ii+2]; /* Yohei */
        } else {
          s[ii+3] = 0.0; /* original */
        }
      }
    }
  }
#endif

  /* calculate normalized hc_limit */
  real nhc_limit = -imin;

  /* calculate height */
  real hm(0.0), hc(0.0), hp(0.0), nhc(0.0);

  for(int ii(imin); ii<=imax; ++ii) {
    hm += stencil[0][ii+3]*gridstencil[0][ii+3];
    hc += stencil[1][ii+3]*gridstencil[1][ii+3];
    hp += stencil[2][ii+3]*gridstencil[2][ii+3];

    nhc += stencil[2][ii+3];
  }

  /* Note: Lopez Eq. 6 correction is not used in 2D */
  if(nhc_limit<nhc && nhc<=(nhc_limit+1.0)) {
    real h_x  = truedir*(hp-hm)/(dm+dp);
    real h_xx = truedir*2.*(dp*hm+dm*hp-(dp+dm)*hc)/hp/hm/(hp+hm);

    kap = mult*h_xx/pow(1.+h_x*h_x,1.5);
  } else {
    flag = 0;
  }

  return;
} 
