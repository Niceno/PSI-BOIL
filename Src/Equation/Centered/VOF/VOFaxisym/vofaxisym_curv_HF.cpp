#include "vofaxisym.h"
#define LOPEZ_COLORING
//#define ONLY_CART
//#define ONLY_CYL
//#define BLEND_CART

/* parameters */
static const int mof(3); /* symmetric stencil is constructed */
static const int nof(1);
static const int majorext(mof*2+1);
static const int minorext(nof*2+1);
/* warning: if stencil size is changed in the minorext direction,
   the kernel must be properly adjusted! */
static const real blending_angle = 37./180.*boil::pi;
static const real n0square = 1.-1./(1.+tan(blending_angle)*tan(blending_angle));
static const real n0 = sqrt(n0square);

static const int iterloop(5); /* 2019.07.09 */
static const real overshoot(0.5);

/******************************************************************************/
void VOFaxisym::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function in 2D axisymmetric.
*     majorext x minorext stencil
*     J.Lopez et al., Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*
*     modification: selection of cells where curvature is calculated is 
*                   modified from Lopez et al.
*
*     modification: adapted for a 2D axisymmetric grid
*
*     in this method, normal vector and bnd values of color are destroyed!
*
*     output: kappa
*******************************************************************************/

  /*-----------------------------------------------------------------+
  |  Step 1: calculate normal vector                                 |
  |  Step 2: define flag=1                                           |
  |  Step 3: calculate curvature if (flag==1 & mof<height<=mof+1)    |
  |  Step 4: calculate curvature near wall                           |
  |  Step 5: extrapolate curvature to iterloop layers                |
  +-----------------------------------------------------------------*/

  /* prepare stencils: one for heights, one for deltas */
  arr2D stencilx, gridstencilx;
  arr2D stencilz, gridstencilz;
  stencilx.resize(minorext);
  gridstencilx.resize(minorext);
  stencilz.resize(minorext);
  gridstencilz.resize(minorext);
  for(auto & s : stencilx) s.resize(majorext);
  for(auto & g : gridstencilx) g.resize(majorext);
  for(auto & s : stencilz) s.resize(majorext);
  for(auto & g : gridstencilz) g.resize(majorext);

  /* Normal vector is used for (i) mMax and (ii) Eq. (2) in J.Lopez et al. */
  if(norm_method_curvature != norm_method_advance) {
    norm(clr,norm_method_curvature,false); /* alpha is not extracted */
  }

  /* detachment treatment = flooding of walls */
  if(detachment_model.initialized()&&detachment_model.detached()) {
    /* this is the only implemented instance atm */
    assert(wall_curv_method==CurvMethod::HFmixedXZ());

    flood(phi,-mult_wall);
    normal_vector_near_bnd(phi,norm_method_curvature);
  }

  /* nx, ny, nz themselves are changed */
  /* mx, my, mz remain unchanged */
  true_norm_vect(nx,ny,nz,nx,ny,nz);

  /* define tempflag */
  tempflag=0;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
#if 0 /* no longer necessary! */
      /* exclude near-wall cells: treated specially */
      if(   (i==si() && iminw) || (i==ei() && imaxw)
         || (k==sk() && kminw) || (k==ek() && kmaxw)
         || dom->ibody().off(i-1,j,k) || dom->ibody().off(i+1,j,k)
         || dom->ibody().off(i,j,k-1) || dom->ibody().off(i,j,k+1)
        ) {
      } else {
#endif
        if((clr[i-1][j][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((clr[i+1][j][k]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((clr[i][j][k-1]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((clr[i][j][k+1]-phisurf)*(clr[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      //}
    }
  }
  tempflag.bnd_update_symmetry(); /* copy on symmetry plane */
  tempflag.exchange();

  /*------------------------+
  |  curvature calculation  |
  +------------------------*/
  kappa=boil::unreal;

  for_ijk(i,j,k) {

    /* calculate curvature only when tempflag=1 */
    if(tempflag[i][j][k]==1) {

      /* select dominant (=majorext) direction of stencil */
      Comp mMax = Comp::undefined();
      /* n points to the liquid */
      real nnx = -nx[i][j][k];
      real nnz = -nz[i][j][k];
      real abs_nx = fabs(nnx);
      real abs_nz = fabs(nnz);
      real max_n;

      if       (abs_nx<n0) {
        mMax = Comp::k();
        max_n = nnz;
      } else if(abs_nz<n0) {
        mMax = Comp::i();
        max_n = nnx;
      }

      if(mMax==Comp::i()) {
        if(!ifull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* fill stencil */
        int imin,imax;
        fill_stencil_x(stencilx,gridstencilx,imin,imax,i,j,k);

        /* calculate heights */
        real hm,hc,hp;
        real nhc; /* for checking if in height range */
        real mult; /* for possible inversion of color */
        calculate_heights(stencilx,gridstencilx,imin,imax,max_n,
                          mult,hm,hc,hp,nhc);
        
        /* calculate curvature, nhc_limit = -imin */
        if(-imin<nhc && nhc<=(-imin+1.0)) {
          /* radius of cylindrical revolution, adjusted for height */
          real xcent_adj = clr.xc(i) + (nhc+imin-0.5)*gridstencilx[nof][mof];

          real kap_cart(0.0), kap_cyl(0.0);
          calculate_curvature_HF_axisymmetric(hm,hc,hp,
                                              clr.dzb(k),clr.dzc(k),clr.dzt(k),
                                              kfull,mult,mMax,xcent_adj,
                                              kap_cart,kap_cyl,
                                              i,j,k);
          kappa[i][j][k] = kap_cart + kap_cyl;
        } else {
          tempflag[i][j][k] = 0;
        }

      } else if(mMax==Comp::k()) {
        if(!kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* fill stencil */
        int kmin,kmax;
        fill_stencil_z(stencilz,gridstencilz,kmin,kmax,i,j,k);

        /* calculate heights */
        real hm,hc,hp;
        real nhc; /* for checking if in height range */
        real mult; /* for possible inversion of color */
        calculate_heights(stencilz,gridstencilz,kmin,kmax,max_n,
                          mult,hm,hc,hp,nhc);

        /* calculate curvature, nhc_limit = -kmin */
        if(-kmin<nhc && nhc<=(-kmin+1.0)) {
          /* radius of cylindrical revolution */
          real xcent = clr.xc(i);

          real kap_cart(0.0), kap_cyl(0.0);
          calculate_curvature_HF_axisymmetric(hm,hc,hp,
                                              clr.dxw(i),clr.dxc(i),clr.dxe(i),
                                              ifull,mult,mMax,xcent,
                                              kap_cart,kap_cyl,
                                              i,j,k);
          kappa[i][j][k] = kap_cart + kap_cyl;
        } else {
          tempflag[i][j][k] = 0;
        }

      } else { /* blending */
        if(!ifull||!kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" while blending; exiting."<<boil::endl;
          exit(0);
        }

        /* fill stencils */
        int imin,imax,kmin,kmax;
        fill_stencil_x(stencilx,gridstencilx,imin,imax,i,j,k);
        fill_stencil_z(stencilz,gridstencilz,kmin,kmax,i,j,k);

        /* calculate heights */
        real hm_x,hc_x,hp_x;
        real nhc_x; /* for checking if in height range */
        real mult_x; /* for possible inversion of color */
        calculate_heights(stencilx,gridstencilx,imin,imax,nnx,
                          mult_x,hm_x,hc_x,hp_x,nhc_x);

        real hm_z,hc_z,hp_z;
        real nhc_z; /* for checking if in height range */
        real mult_z; /* for possible inversion of color */
        calculate_heights(stencilz,gridstencilz,kmin,imax,nnz,
                          mult_z,hm_z,hc_z,hp_z,nhc_z);

        bool flag_x = -imin<nhc_x && nhc_x<=(-imin+1.0);
        bool flag_z = -kmin<nhc_z && nhc_z<=(-kmin+1.0);
        bool flag_x_adj = -imin-overshoot<nhc_x && nhc_x<=(-imin+1.0+overshoot);
        bool flag_z_adj = -kmin-overshoot<nhc_z && nhc_z<=(-kmin+1.0+overshoot);
#if 0
        bool flag = flag_x & flag_z;
        boil::oout<<nhc_x<<" "<<nhc_z<<" | "<<nhc_x+imin<<" "<<nhc_z+kmin<<" | "
                  <<flag_x<<" "<<flag_z<<" "<<flag
                  <<boil::endl;
#endif
        real kap_x_cyl(0.0), kap_z_cyl(0.0);
        real kap_x_cart(0.0), kap_z_cart(0.0);
        
        if(!flag_x&&!flag_z) {
          tempflag[i][j][k] = 0;
        } else {
          /* blending factor for z-component */
          real bfactor = ( flag_x_adj& flag_z_adj)*(abs_nz*abs_nz - n0square)/(1.-2.*n0square)
                       + ( flag_x_adj&!flag_z_adj)*0.0
                       + (!flag_x_adj& flag_z_adj)*1.0;
          if(flag_x_adj) {
            /* radius of cylindrical revolution, adjusted for height */
            real xcent_x = clr.xc(i) + (nhc_x+imin-0.5)*gridstencilx[nof][mof];

            calculate_curvature_HF_axisymmetric(hm_x,hc_x,hp_x,
                                                clr.dzb(k),clr.dzc(k),clr.dzt(k),
                                                kfull,mult_x,Comp::i(),xcent_x,
                                                kap_x_cart,kap_x_cyl,
                                                i,j,k);
            kap_x_cyl *= 1.0-bfactor;
#ifdef BLEND_CART
            kap_x_cart *= 1.0-bfactor;
#endif
          }
          if(flag_z_adj) {
            /* radius of cylindrical revolution */
            real xcent_z = clr.xc(i);

            calculate_curvature_HF_axisymmetric(hm_z,hc_z,hp_z,
                                                clr.dxw(i),clr.dxc(i),clr.dxe(i),
                                                ifull,mult_z,Comp::k(),xcent_z,
                                                kap_z_cart,kap_z_cyl,
                                                i,j,k);
            kap_z_cyl *= bfactor;
#ifdef BLEND_CART
            kap_z_cart *= bfactor;
#endif
          }

          kappa[i][j][k] = kap_x_cyl + kap_z_cyl;
#ifndef BLEND_CART
          /* majorext direction still used for cartesian curvature */
          if(abs_nx<abs_nz) {
            kappa[i][j][k] += kap_z_cart;
          } else {
            kappa[i][j][k] += kap_x_cart;
          }
#else
          kappa[i][j][k] += kap_x_cart + kap_z_cart;
#endif

        }
      } /* blending */

      //boil::oout<<i<<" "<<k<<" "<<mMax<<" "<<kappa[i][j][k]<<boil::endl;

    } /* tempflag = 1 */
  } /* for ijk */

  kappa.bnd_update_symmetry(); /* copy on symmetry plane */
  tempflag.bnd_update_symmetry(); /* copy on symmetry plane */
  kappa.exchange();
  tempflag.exchange();

  /* near-wall calculations */
  if       (wall_curv_method==CurvMethod::DivNorm()) {
    boil::oout<<"VOFaxisym::curv_HF: Obsolete, underdeveloped code. Exiting."
              <<boil::endl;
    exit(0);
    //bdcurv();
  } else if(wall_curv_method==CurvMethod::HFmixedXZ()) {
    insert_bc_curv_HFmixed(clr,Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::none()) {
  } else {
    /* default */
    boil::oout<<"VOF::curvHF: Wall curvature calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }

  /* extrapolate kappa */
  stmp  = kappa;
  tempflag2 = tempflag;
  
  for(int iloop=1; iloop<iterloop; iloop++) { 
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(tempflag[i][j][k]==0) {
        /* at this point, the near-wall tempflag must be properly 0 or 1 */
        int inb = std::min(1,tempflag[i-1][j][k]) + std::min(1,tempflag[i+1][j][k])
                + std::min(1,tempflag[i][j][k-1]) + std::min(1,tempflag[i][j][k+1]);
        if(inb >= 1) {
          stmp[i][j][k] = (real(std::min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                        +  real(std::min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                        +  real(std::min(1,tempflag[i][j][k-1])) * kappa[i][j][k-1]
                        +  real(std::min(1,tempflag[i][j][k+1])) * kappa[i][j][k+1])
                        /real(inb);
          tempflag2[i][j][k] = 2;  /* tempflag=2 for extrapolated */
        }
      }
    }
    stmp.bnd_update_symmetry(); /* copy on symmetry plane */
    tempflag2.bnd_update_symmetry(); /* copy on symmetry plane */
    stmp.exchange();
    tempflag2.exchange();
    kappa = stmp;
    tempflag = tempflag2;
  }

#if 0
  //if(time->current_step()==1) {
  //if(time->current_step()%100==0) {
  /* visualize tempflag */
  boil::plot->plot(clr,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=tempflag[i][j][k];
  }
  boil::plot->plot(clr,kappa,stmp, "clr-kappa-tempflag", time->current_step());
  exit(0);
  //}
#endif
  
  return;
}

/*----------------+
|  fill stencils  |
+----------------*/
void VOFaxisym::fill_stencil_x(arr2D & stencil, arr2D & gridstencil,
                               int & imin, int & imax,
                               const int i, const int j, const int k) {

  /* check stencil size */
  imin=-mof;
  imax=mof;  /* normal stencil size 2*mof+1 */

  /* limit stencil size for cut-stencil */
  if(iminc) imin=std::max(-mof,si()-i);
  if(imaxc) imax=std::min( mof,ei()-i);

  for(int ii(2); ii<=mof; ++ii) {
    if(dom->ibody().off(i-ii,j,k)) {
      imin=-ii+1;
      break;
    }
  }
  for(int ii(2); ii<=mof; ++ii) {
    if(dom->ibody().off(i+ii,j,k)) {
      imax=ii-1;
      break;
    }
  }

  /* fill stencil */
  for(int kk(-nof); kk<=nof; ++kk) {
    /* stencil is reset with negative values representing undefined vals */
    for(int ii(-mof); ii<=mof; ++ii) {
      gridstencil[kk+nof][ii+mof] = 0.0;
      stencil[kk+nof][ii+mof] = -1.0;
    }
    for(int ii(imin); ii<=imax; ++ii) {
      gridstencil[kk+nof][ii+mof] = clr.dxc(i+ii);
      stencil[kk+nof][ii+mof] = std::min(1.0,std::max(0.0,clr[i+ii][j][k+kk]));
    }
  }

  return;
}

void VOFaxisym::fill_stencil_z(arr2D & stencil, arr2D & gridstencil,
                               int & kmin, int & kmax,
                               const int i, const int j, const int k) {
  /* check stencil size */
  kmin=-mof;
  kmax=mof;  /* normal stencil size 2*mof+1 */

  /* limit stencil size for cut-stencil */
  if(kminc) kmin=std::max(-mof,sk()-k);
  if(kmaxc) kmax=std::min( mof,ek()-k);

  for(int kk(2); kk<=mof; ++kk) {
    if(dom->ibody().off(i,j,k-kk)) {
      kmin=-kk+1;
      break;
    }
  }
  for(int kk(2); kk<=mof; ++kk) {
    if(dom->ibody().off(i,j,k+kk)) {
      kmax=kk-1;
      break;
    }
  }

  /* fill stencil */
  for(int ii(-nof); ii<=nof; ++ii) {
    /* stencil is reset with negative values representing undefined vals */
    for(int kk(-mof); kk<=mof; ++kk) {
      gridstencil[ii+nof][kk+mof] = 0.0;
      stencil[ii+nof][kk+mof] = -1.0;
    }
    for(int kk(kmin); kk<=kmax; ++kk) {
      gridstencil[ii+nof][kk+mof] = clr.dzc(k+kk);
      stencil[ii+nof][kk+mof] = std::min(1.0,std::max(0.0,clr[i+ii][j][k+kk]));
    }
  }

  return;
}

/*---------------------+
|  height calculation  |
+---------------------*/
void VOFaxisym::calculate_heights(arr2D & stencil, const arr2D & gridstencil,
                                  const int imin, const int imax,
                                  const real max_n, real & mult,
                                  real & hm, real & hc, real & hp, 
                                  real & nhc) {

  /* color is inverted if normal vector points in negative dir */
  mult = 1.0;
  if(max_n<0.0) {
    mult = -1.0;
    for(auto & s : stencil) {
      for(int ii(imin); ii<=imax; ++ii) {
        s[ii+mof] = 1.0-s[ii+mof];
      }
    }
  }

  /* Calculate Eq. (2) in J.Lopez et al. */
  /* Warning! if max_n is low accuracy (wrong sign), this algorithm
     doesn't work !!!  (Yohei) */
#ifdef LOPEZ_COLORING
  for(auto & s : stencil) {
    for(int ii(-1); ii>=imin; --ii) {
      if( (s[ii+mof]-s[ii+mof+1]) < 0.0) {
        if(  (s[ii+mof]-phisurf)*(s[ii+mof+1]-phisurf) > 0.0
           &&(s[ii+mof]-phisurf)*(s[mof     ]-phisurf) > 0.0) {
          s[ii+mof] = s[ii+mof+1]; /* Yohei */
        } else {
          s[ii+mof] = 1.0; /* original */
        }
      }
    }

    for(int ii( 1); ii<=imax; ++ii) {
      if( (s[ii+mof]-s[ii+mof-1]) > 0.0) {
        if(  (s[ii+mof]-phisurf)*(s[ii+mof-1]-phisurf) > 0.0
           &&(s[ii+mof]-phisurf)*(s[mof     ]-phisurf) > 0.0) {
          s[ii+mof] = s[ii+mof-1]; /* Yohei */
        } else {
          s[ii+mof] = 0.0; /* original */
        }
      }

    }
  }
#endif

  hm = hc = hp = nhc = 0.0;
  for(int ii(imin); ii<=imax; ++ii) {
    hm += stencil[0][ii+mof]*gridstencil[0][ii+mof];
    hc += stencil[1][ii+mof]*gridstencil[1][ii+mof];
    hp += stencil[2][ii+mof]*gridstencil[2][ii+mof];

    nhc += stencil[nof][ii+mof];
  }
  
  return;
}

/*-----------------------+
|  curvature calculation |
+-----------------------*/
void VOFaxisym::calculate_curvature_HF_axisymmetric(
                               const real hm, const real hc, const real hp,
                               const real dm, const real dc, const real dp,
                               const bool truedir, const real mult,
                               const Comp mcomp, const real xcent,
                               real & kap_cart, real & kap_cyl,
                               const int i, const int j, const int k) {

/* Note: Lopez Eq. (6) correction is not used in 2D */
  real h_1c = truedir*(hp-hm)/(dm+dp);
  real h_1u = truedir*(hp-hc)/dp;
  real h_1d = truedir*(hc-hm)/dm;
  real h_11 = truedir*2.*(dp*hm+dm*hp-(dp+dm)*hc)/dp/dm/(dp+dm);

  /* Cartesian contribution */
  kap_cart = h_11/pow(1.+h_1c*h_1c,1.5);

  /* cylindrical contribution -> depends on stencil orientation */
  if(mcomp==Comp::k()) {
#if 0
    kap_cyl  = 1.0/xcent*h_1c/sqrt(1.+h_1c*h_1c);
#else
    kap_cyl  = h_1u/sqrt(1.+h_1u*h_1u);
    kap_cyl += h_1d/sqrt(1.+h_1d*h_1d);

    kap_cyl /= 2.*xcent;
#endif
  } else {
#if 1
    kap_cyl = -1./xcent * 1./sqrt(1.+h_1c*h_1c);
#else /* this cannot reproduce zero curvature at inflexion */
    kap_cyl  = -1./sqrt(1.+h_1u*h_1u);
    kap_cyl += -1./sqrt(1.+h_1d*h_1d);

    kap_cyl /= 2.*xcent;
#endif
  }

  /* under this convention, bubbles have negative curvature */
  kap_cart *= -mult;
  kap_cyl *= -mult;

#ifdef ONLY_CART
  kap_cyl *= 0.0;
#elif defined ONLY_CYL
  kap_cart *= 0.0;
#endif

  //boil::oout<<i<<" "<<k<<" | "<<h_1d<<" "<<h_1c<<" "<<h_1u<<" "<<h_11<<" | "<<mcomp<<" | "<<kap_cart<<" "<<kap_cyl<<" "<<kap_cart+kap_cyl<<boil::endl;

  return;
} 
