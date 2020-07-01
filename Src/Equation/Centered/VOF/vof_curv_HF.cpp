#include "vof.h"
#define LOPEZ_COLORING
#define LOPEZ_SMOOTHING

/******************************************************************************/
void VOF::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function in 2D axisymmetric.
*     majorext x minorext x minorext stencil
*     J.Lopez et al., Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*
*     modification: selection of cells where curvature is calculated is 
*                   modified from Lopez et al.
*
*     in this method, normal vector and bnd values of color are destroyed!
*
*     output: kappa
*******************************************************************************/

  const int & mof = hf_set.mof;
  const int & nof = hf_set.nof;
  const int & majorext = hf_set.majorext;
  const int & minorext = hf_set.minorext;

  /*-----------------------------------------------------------------+
  |  Step 1: calculate normal vector                                 |
  |  Step 2: define flag=1                                           |
  |  Step 3: calculate curvature if (flag==1 & mof<height<=mof+1)    |
  |  Step 4: calculate curvature near wall                           |
  |  Step 5: extrapolate curvature to iterloop layers                |
  +-----------------------------------------------------------------*/

  /* prepare stencils: one for heights, one for dx */
  arr3D stencil, gridstencil;
  stencil.resize(minorext);
  gridstencil.resize(minorext);
  for(auto & t : stencil) {
    t.resize(minorext);
    for(auto & s : t) {
      s.resize(majorext);
    }
  }
  for(auto & t : gridstencil) {
    t.resize(minorext);
    for(auto & g : t) {
      g.resize(majorext);
    }
  }

  /* Normal vector is used for (i) mMax and (ii) Eq. (2) in J.Lopez et al. */
  if(norm_method_curvature != norm_method_advance) {
    norm(phi,norm_method_curvature,false); /* alpha is not extracted */
  }

  /* detachment treatment = flooding of walls */
  if(detachment_model.initialized()&&detachment_model.detached()) {
    /* this is the only implemented instance atm */
    if(wall_curv_method==CurvMethod::HFparallelXZ()) {
      flood(phi,-mult_wall);
      normal_vector_near_bnd(phi,norm_method_curvature);
    }
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
         || (j==sj() && jminw) || (j==ej() && jmaxw)
         || (k==sk() && kminw) || (k==ek() && kmaxw)
         || dom->ibody().off(i-1,j,k) || dom->ibody().off(i+1,j,k)
         || dom->ibody().off(i,j-1,k) || dom->ibody().off(i,j+1,k)
         || dom->ibody().off(i,j,k-1) || dom->ibody().off(i,j,k+1)
         || (i==si()+1 && iminw) || (i==ei()-1 && imaxw)
         || (j==sj()+1 && jminw) || (j==ej()-1 && jmaxw)
         || (k==sk()+1 && kminw) || (k==ek()-1 && kmaxw)
         || dom->ibody().off(i-2,j,k) || dom->ibody().off(i+2,j,k)
         || dom->ibody().off(i,j-2,k) || dom->ibody().off(i,j+2,k)
         || dom->ibody().off(i,j,k-2) || dom->ibody().off(i,j,k+2)
        ) {
      } else {
#endif
        if((phi[i-1][j][k]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((phi[i+1][j][k]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((phi[i][j-1][k]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((phi[i][j+1][k]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((phi[i][j][k-1]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
        if((phi[i][j][k+1]-phisurf)*(phi[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
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
      Comp mMax;
      /* n points to the liquid */
      real nnx = -nx[i][j][k];
      real nny = -ny[i][j][k];
      real nnz = -nz[i][j][k];
      real abs_nx = fabs(nnx);
      real abs_ny = fabs(nny);
      real abs_nz = fabs(nnz);
      real max_n;

      if(abs_nx<abs_ny) {
        if(abs_ny<abs_nz) {
          mMax = Comp::k();
          max_n = nnz;
        } else {
          mMax = Comp::j();
          max_n = nny;
        }
      } else {
        if(abs_nx<abs_nz) {
          mMax = Comp::k();
          max_n = nnz;
        } else {
          mMax = Comp::i();
          max_n = nnx;
        }
      }

      if(mMax==Comp::i()) {
        if(!ifull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int imin=-mof;
        int imax= mof;  /* normal stencil size 2*mof+1 */

        /* limit stencil size for cut-stencil */
        if(iminc) imin=std::max(-mof,si()-i);
        if(imaxc) imax=std::min( mof,ei()-i);

        /* because of ghost grid in solid/walls */
        int imin_grid = imin;
        int imax_grid = imax;
        
        /* update at walls */
        if( iminw && -mof<si()-i ) {
          imin = si()-i-1;
          imin_grid = imin+1;
        }
        if( imaxw &&  mof>ei()-i ) {
          imax = ei()-i+1;
          imax_grid = imax-1;
        }

        for(int ii(1); ii<=mof; ++ii) {
          if(dom->ibody().off(i-ii,j,k)) {
            imin = -ii;
            imin_grid = imin+1;
            break;
          }
        }
        for(int ii(1); ii<=mof; ++ii) {
          if(dom->ibody().off(i+ii,j,k)) {
            imax = ii;
            imax_grid = imax-1;
            break;
          }
        }

        /* fill stencil */
        for(int jj(-nof); jj<=nof; ++jj) {
          for(int kk(-nof); kk<=nof; ++kk) {
            /* stencil is reset with negative values representing undefined vals */
            for(int ii(-mof); ii<=mof; ++ii) {
              gridstencil[jj+nof][kk+nof][ii+mof] = 0.0;
              stencil[jj+nof][kk+nof][ii+mof] = -1.0;
            }
            for(int ii(imin); ii<=imax; ++ii) {
              gridstencil[jj+nof][kk+nof][ii+mof] = phi.dxc(i+ii);
              stencil[jj+nof][kk+nof][ii+mof] = std::min(1.0,std::max(0.0,phi[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(imin_grid != imin) {
          for(int jj(-nof); jj<=nof; ++jj) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[jj+nof][kk+nof][imin+mof] = phi.dxc(i+imin_grid);
            }
          }
        }
        if(imax_grid != imax) {
          for(int jj(-nof); jj<=nof; ++jj) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[jj+nof][kk+nof][imax+mof] = phi.dxc(i+imax_grid);
            }
          }
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,
                       imin,imax,
                       phi.dys(j),phi.dyc(j),phi.dyn(j),
                       phi.dzb(k),phi.dzc(k),phi.dzt(k),
                       max_n,jfull,kfull,
                       kappa[i][j][k],tempflag[i][j][k]);
                       //i,j,k);

      } else if(mMax==Comp::j()) {
        if(!jfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int jmin=-mof;
        int jmax= mof;  /* normal stencil size 2*mof+1 */

        /* limit stencil size for cut-stencil */
        if(jminc) jmin=std::max(-mof,sj()-j);
        if(jmaxc) jmax=std::min( mof,ej()-j);

        /* because of ghost grid in solid/walls */
        int jmin_grid = jmin;
        int jmax_grid = jmax;
        
        /* update at walls */
        if( jminw && -mof<sj()-j ) {
          jmin = sj()-j-1;
          jmin_grid = jmin+1;
        }
        if( jmaxw &&  mof>ej()-j ) {
          jmax = ej()-j+1;
          jmax_grid = jmax-1;
        }

        for(int jj(1); jj<=mof; ++jj) {
          if(dom->ibody().off(i,j-jj,k)) {
            jmin = -jj;
            jmin_grid = jmin+1;
            break;
          }
        }
        for(int jj(1); jj<=mof; ++jj) {
          if(dom->ibody().off(i,j+jj,k)) {
            jmax = jj;
            jmax_grid = jmax-1;
            break;
          }
        }

        /* fill stencil */
        for(int ii(-nof); ii<=nof; ++ii) {
          for(int kk(-nof); kk<=nof; ++kk) {
            /* stencil is reset with negative values representing undefined vals */
            for(int jj(-mof); jj<=mof; ++jj) {
              gridstencil[ii+nof][kk+nof][jj+mof] = 0.0;
              stencil[ii+nof][kk+nof][jj+mof] = -1.0;
            }
            for(int jj(jmin); jj<=jmax; ++jj) {
              gridstencil[ii+nof][kk+nof][jj+mof] = phi.dyc(j+jj);
              stencil[ii+nof][kk+nof][jj+mof] = std::min(1.0,std::max(0.0,phi[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(jmin_grid != jmin) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[ii+nof][kk+nof][jmin+mof] = phi.dyc(j+jmin_grid);
            }
          }
        }
        if(jmax_grid != jmax) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[ii+nof][kk+nof][jmax+mof] = phi.dyc(j+jmax_grid);
            }
          }
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,
                       jmin,jmax,
                       phi.dxw(i),phi.dxc(i),phi.dxe(i),
                       phi.dzb(k),phi.dzc(k),phi.dzt(k),
                       max_n,ifull,kfull,
                       kappa[i][j][k],tempflag[i][j][k]);
                       //i,j,k);

      } else { 
        if(!kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int kmin=-mof;
        int kmax= mof;  /* normal stencil size 2*mof+1 */ 

        /* limit stencil size for cut-stencil */
        if(kminc) kmin=std::max(-mof,sk()-k);
        if(kmaxc) kmax=std::min( mof,ek()-k);

        /* because of ghost grid in solid/walls */
        int kmin_grid = kmin;
        int kmax_grid = kmax;
        
        /* update at walls */
        if( kminw && -mof<sk()-k ) {
          kmin = sk()-k-1;
          kmin_grid = kmin+1;
        }
        if( kmaxw &&  mof>ek()-k ) {
          kmax = ek()-k+1;
          kmax_grid = kmax-1;
        }

        for(int kk(1); kk<=mof; ++kk) {
          if(dom->ibody().off(i,j,k-kk)) {
            kmin = -kk;
            kmin_grid = kmin+1;
            break;
          }
        }
        for(int kk(1); kk<=mof; ++kk) {
          if(dom->ibody().off(i,j,k+kk)) {
            kmax = kk;
            kmax_grid = kmax-1;
            break;
          }
        }

        /* fill stencil */
        for(int ii(-nof); ii<=nof; ++ii) {
          for(int jj(-nof); jj<=nof; ++jj) {
            /* stencil is reset with negative values representing undefined vals */
            for(int kk(-mof); kk<=mof; ++kk) {
              gridstencil[ii+nof][jj+nof][kk+mof] = 0.0;
              stencil[ii+nof][jj+nof][kk+mof] = -1.0;
            }
            for(int kk(kmin); kk<=kmax; ++kk) {
              gridstencil[ii+nof][jj+nof][kk+mof] = phi.dzc(k+kk);
              stencil[ii+nof][jj+nof][kk+mof] = std::min(1.0,std::max(0.0,phi[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(kmin_grid != kmin) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int jj(-nof); jj<=nof; ++jj) {
              gridstencil[ii+nof][jj+nof][kmin+mof] = phi.dzc(k+kmin_grid);
            }
          }
        }
        if(kmax_grid != kmax) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int jj(-nof); jj<=nof; ++jj) {
              gridstencil[ii+nof][jj+nof][kmax+mof] = phi.dzc(k+kmax_grid);
            }
          }
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,
                       kmin,kmax,
                       phi.dxw(i),phi.dxc(i),phi.dxe(i),
                       phi.dys(j),phi.dyc(j),phi.dyn(j),
                       max_n,ifull,jfull,
                       kappa[i][j][k],tempflag[i][j][k]);
                       //i,j,k);

      }

    } /* tempflag = 1 */
  } /* for ijk */

  kappa.bnd_update_symmetry(); /* copy on symmetry plane */
  tempflag.bnd_update_symmetry(); /* copy on symmetry plane */
  kappa.exchange();
  tempflag.exchange();

  /* near-wall calculations */
  if       (wall_curv_method==CurvMethod::DivNorm()) {
    insert_bc_curv_divnorm();
  } else if(wall_curv_method==CurvMethod::HFmixedXZ()) {
    insert_bc_curv_HFmixed(phi,Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::HFparallelXZ()) {
    insert_bc_curv_HFparallel(phi,Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::HFnormalXZ()) {
    insert_bc_curv_HFnormal(phi,Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::HFnormalZ()) {
    insert_bc_curv_HFnormal(phi,Comp::k(),Sign::neg());
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
  
  //if(false) {
  for(int iloop=1; iloop<hf_set.iterloop; iloop++) { /* 2019.07.09 */
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) continue;
      if(tempflag[i][j][k]==0) {
        /* at this point, the near-wall tempflag must be properly 0 or 1 */
        int inb = std::min(1,tempflag[i-1][j][k]) + std::min(1,tempflag[i+1][j][k])
                + std::min(1,tempflag[i][j-1][k]) + std::min(1,tempflag[i][j+1][k])
                + std::min(1,tempflag[i][j][k-1]) + std::min(1,tempflag[i][j][k+1]);
        if(inb >= 1) {
          stmp[i][j][k] = (real(std::min(1,tempflag[i-1][j][k])) * kappa[i-1][j][k]
                        +  real(std::min(1,tempflag[i+1][j][k])) * kappa[i+1][j][k]
                        +  real(std::min(1,tempflag[i][j-1][k])) * kappa[i][j-1][k]
                        +  real(std::min(1,tempflag[i][j+1][k])) * kappa[i][j+1][k]
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
  if(time->current_step()==1) {
  //if(time->current_step()%100==0) {
  /* visualize tempflag */
  boil::plot->plot(phi,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  for_ijk(i,j,k){
    stmp[i][j][k]=tempflag[i][j][k];
  }
  boil::plot->plot(phi,kappa,stmp, "clr-kappa-tempflag", time->current_step());
  exit(0);
  } 
#endif

  return;
}

/*-------------------------------+
|  curvature calculation kernel  |
+-------------------------------*/
/* inputs:
   - stencil of color values
   - stencil of grid spacings in the dominant direction
   - starting and ending indices in the dominant direction
   - grid spacings in the first minorext direction
   - grid spacings in the second minorext direction
   - normal vector component in the dominant direction, pointing to gas
   - bool indicator if first minorext direction is a non-pseudo direction
   - bool indicator if second minorext direction is a non-pseudo direction
   
   outputs:
   - kappa value (if computed)
   - tempflag value (if kappa not computed)
*/
void VOF::curv_HF_kernel(arr3D & stencil, const arr3D & gridstencil,
                         const int imin, const int imax,
                         const real d1m, const real d1c, const real d1p,
                         const real d2m, const real d2c, const real d2p,
                         const real max_n,
                         const bool truedir1, const bool truedir2,
                         real & kap, int & flag) {
                         //const int i, const int j, const int k) {

  const int & mof = hf_set.mof;
  const int & nof = hf_set.nof;

  /* color is inverted if normal vector points in negative dir */
  real mult(1.0);
  if(max_n<0.0) {
    mult = -1.0;
    for(auto & t : stencil) {
      for(auto & s : t) {
        for(int ii(imin); ii<=imax; ++ii) {
          s[ii+mof] = 1.0-s[ii+mof];
        }
      }
    }
  }

  /* Calculate Eq. (2) in J.Lopez et al. */
  /* Warning! if max_n is low accuracy (wrong sign), this algorithm
     doesn't work !!!  (Yohei) */
#ifdef LOPEZ_COLORING
  for(auto & t : stencil) {
    for(auto & s : t) {
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
  }
#endif

  /* calculate normalized hc_limit */
  real nhc_limit = -imin;

  /* calculate height */
  real hmm(0.0), hcm(0.0), hpm(0.0);
  real hmc(0.0), hcc(0.0), hpc(0.0), nhc(0.0);
  real hmp(0.0), hcp(0.0), hpp(0.0);

  for(int ii(imin); ii<=imax; ++ii) {
    hmm += stencil[0][0][ii+mof]*gridstencil[0][0][ii+mof]; 
    hcm += stencil[1][0][ii+mof]*gridstencil[1][0][ii+mof]; 
    hpm += stencil[2][0][ii+mof]*gridstencil[2][0][ii+mof]; 

    hmc += stencil[0][1][ii+mof]*gridstencil[0][1][ii+mof]; 
    hcc += stencil[1][1][ii+mof]*gridstencil[1][1][ii+mof]; 
    hpc += stencil[2][1][ii+mof]*gridstencil[2][1][ii+mof]; 

    hmp += stencil[0][2][ii+mof]*gridstencil[0][2][ii+mof]; 
    hcp += stencil[1][2][ii+mof]*gridstencil[1][2][ii+mof]; 
    hpp += stencil[2][2][ii+mof]*gridstencil[2][2][ii+mof]; 

    nhc += stencil[nof][nof][ii+mof];
  }

  if(nhc_limit<nhc && nhc<=(nhc_limit+1.0)) {

    kap = calculate_curvature_HF(hmm, hmc, hmp,
                                 hcm, hcc, hcp,
                                 hpm, hpc, hpp,
                                 d1m, d1c, d1p,
                                 d2m, d2c, d2p,
                                 truedir1, truedir2,
                                 mult, max_n); 

#if 0
    boil::oout<<i<<" "<<j<<" "<<k<<" | "
              <<hmm<<" "<<hcm<<" "<<hpm<<" | "
              <<hmc<<" "<<hcc<<" "<<hpc<<" "
              <<hmp<<" "<<hcp<<" "<<hpp<<" "
              //<<h_1<<" "<<h_2<<" "<<h_11<<" "<<h_22<<" "<<h_12
              <<boil::endl;
#endif
  } else {
    flag = 0;
  }

  return;
} 

real VOF::calculate_curvature_HF(
                               const real hmm, const real hmc, const real hmp,
                               const real hcm, const real hcc, const real hcp,
                               const real hpm, const real hpc, const real hpp,
                               const real d1m, const real d1c, const real d1p,
                               const real d2m, const real d2c, const real d2p,
                               const bool truedir1, const bool truedir2,
                               const real mult, const real max_n) {
  const real & theta_crit = hf_set.theta_crit;
  real theta = acos(mult*max_n);
  real g = 0.0;
  if(theta < theta_crit) { /* Eq. (6) Lopez (2009) Comp Methods Appl */
    g = 0.2;
  }

#ifndef LOPEZ_SMOOTHING
  real h_1  = truedir1*(hpc-hmc)/(d1m+d1p);
  real h_11 = truedir1*2.*(d1p*hmc+d1m*hpc-(d1p+d1m)*hcc)/d1p/d1m/(d1p+d1m);
  real h_2  = truedir2*(hcp-hcm)/(d2m+d2p);
  real h_22 = truedir2*2.*(d2p*hcm+d2m*hcp-(d2p+d2m)*hcc)/d2p/d2m/(d2p+d2m);
#else /* Lopez correction */
  real h_1  = truedir1*(g*(hpm-hmm)+(hpc-hmc)+g*(hpp-hmp))/(d1m+d1p)/(1.+2.*g);
  real h_11 = truedir1*4.*(
                g*(hmm+hpm-2.*hcm)
               +  (hmc+hpc-2.*hcc)
               +g*(hmp+hpp-2.*hcp)
              )/(d1p+d1m)/(d1p+d1m)/(1.+2.*g);
    
  real h_2  = truedir2*(g*(hmp-hmm)+(hcp-hcm)+g*(hpp-hpm))/(d2m+d2p)/(1.+2.*g);
  real h_22 = truedir2*2.*(
                g*(d2p*hpm+d2m*hpp-(d2p+d2m)*hpc)
               +  (d2p*hcm+d2m*hcp-(d2p+d2m)*hcc)
               +g*(d2p*hmm+d2m*hmp-(d2p+d2m)*hmc)
              )/d2p/d2m/(d2p+d2m)/(1.+2.*g);
#endif
  real h_12 = truedir1*truedir2*(hpp-hpm-hmp+hmm)/(d1p+d1m)/(d2p+d2m);
    
  /* under this convention, bubbles have negative curvature */
  return -mult
         * (h_11 + h_22 + h_11*h_2*h_2 + h_22*h_1*h_1 - 2.0*h_12*h_1*h_2)
         / pow(1.+h_1*h_1 + h_2*h_2,1.5);
} 
