#include "vof.h"

/******************************************************************************/
void VOF::curv_HF() {
/***************************************************************************//**
*  \brief Calculate curvature using height function in 3D cartesian.
*     majorext x minorext x minorext stencil
*     J.Lopez et al., Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*
*     modification: selection of cells where curvature is calculated is 
*                   modified from Lopez et al.
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

  /* ibody/wall */
  arr2D wall_indicator;
  wall_indicator.resize(minorext);
  for(auto & wi : wall_indicator) {
    wi.resize(minorext);
  }

  /* Normal vector is used for (i) mMax and (ii) Eq. (2) in J.Lopez et al. */
  if(norm_method_curvature != norm_method_advance) {
    norm(color(),norm_method_curvature,false); /* alpha is not extracted */
  }

#if 0
  /* detachment treatment = flooding of walls */
  if(detachment_model.initialized()&&detachment_model.detached()) {
    /* this is the only implemented instance atm */
    if(wall_curv_method==CurvMethod::HFparallelXZ()) {
      flood(color(),-mult_wall);
      normal_vector_near_bnd(color(),norm_method_curvature);
    }
  }
#endif

  /* nx, ny, nz themselves are changed */
  /* mx, my, mz remain unchanged */
  true_norm_vect(nx,ny,nz,nx,ny,nz);

  /* define tempflag */
  tempflag=0;
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      if((color()[i-1][j][k]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      if((color()[i+1][j][k]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      if((color()[i][j-1][k]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      if((color()[i][j+1][k]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      if((color()[i][j][k-1]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
      if((color()[i][j][k+1]-phisurf)*(color()[i][j][k]-phisurf)<=0.0){ tempflag[i][j][k]=1;}
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

      /* near wall, normal vector might be of low quality: enforced dir */
      std::array<bool,3> nearwall = {false,false,false};
      if(  (i==si()&&bflag_struct.iminw)
         ||(i==ei()&&bflag_struct.imaxw)
         ||dom->ibody().off(i-1,j,k)
         ||dom->ibody().off(i+1,j,k)) {
        nearwall[0] = true;
      }
      if(  (j==sj()&&bflag_struct.jminw)
         ||(j==ej()&&bflag_struct.jmaxw)
         ||dom->ibody().off(i,j-1,k)
         ||dom->ibody().off(i,j+1,k)) {
        nearwall[1] = true;
      }
      if(  (k==sk()&&bflag_struct.kminw)
         ||(k==ek()&&bflag_struct.kmaxw)
         ||dom->ibody().off(i,j,k-1)
         ||dom->ibody().off(i,j,k+1)) {
        nearwall[2] = true;
      }

      /* next to at least one wall */
      if(nearwall[0]||nearwall[1]||nearwall[2]) {
        if(wall_curv_dir!=Comp::undefined()) {
          mMax = wall_curv_dir;
        /* only apply adaptive method when next to one wall */
        } else if( (nearwall[0]+nearwall[1]+nearwall[2])<2 ) {
          /* local cangle-based selection */
          real cng_c = cangle(i,j,k);
          if(cng_c>boil::pi/4.&&cng_c<3.*boil::pi/4.) {
            /* i dont think anything special has to be done in this case - LB */
          } else {
            /* choose according to which wall it is */
            mMax = nearwall[0] ? Comp::i()
                 : nearwall[1] ? Comp::j()
                               : Comp::k();
          }
        }

        if(mMax==Comp::i()) {
          max_n = nnx;
        } else if(mMax==Comp::j()) {
          max_n = nny;
        } else {
          max_n = nnz;
        }
      }

      /********************* height construction */
      if(mMax==Comp::i()) {
        if(!bflag_struct.ifull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int imin=-mof;
        int imax= mof;  /* normal stencil size 2*mof+1 */

        /* limit stencil size for cut-stencil */
        if(bflag_struct.iminc) imin=std::max(-mof,si()-i);
        if(bflag_struct.imaxc) imax=std::min( mof,ei()-i);

        /* because of ghost grid in solid/walls */
        int imin_grid = imin;
        int imax_grid = imax;
        
        /* update at walls */
        if( bflag_struct.iminw && -mof<si()-i ) {
          imin = si()-i-1;
          imin_grid = imin+1;
        }
        if( bflag_struct.imaxw &&  mof>ei()-i ) {
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

        /* cangle values */
        real cang_m(-10.), cang_p(-10.);

        if(imax_grid!=imax) {
          cang_p = cangle(i+imax_grid,j,k);
        }
        if(imin_grid!=imin) {
          cang_m = cangle(i+imin_grid,j,k);
        }

        /* fill wall_indicator */
        for(int jj(-nof); jj<=nof; ++jj) {
          for(int kk(-nof); kk<=nof; ++kk) {
            /* central position has the contact angle */
            if(jj==0&&kk==0) {
              wall_indicator[jj+nof][kk+nof] = cangle(i,j,k); 
            /* indicator is reset */
            } else if(dom->ibody().off(i,j+jj,k+kk)) {
              wall_indicator[jj+nof][kk+nof] = 0.0;
            } else {
              wall_indicator[jj+nof][kk+nof] = 1.0;
            }
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
              gridstencil[jj+nof][kk+nof][ii+mof] = color().dxc(i+ii);
              stencil[jj+nof][kk+nof][ii+mof] = std::min(1.0,std::max(0.0,color()[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(imin_grid != imin) {
          for(int jj(-nof); jj<=nof; ++jj) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[jj+nof][kk+nof][imin+mof] = color().dxc(i+imin_grid);
            }
          }
        }
        if(imax_grid != imax) {
          for(int jj(-nof); jj<=nof; ++jj) {
            for(int kk(-nof); kk<=nof; ++kk) {
              gridstencil[jj+nof][kk+nof][imax+mof] = color().dxc(i+imax_grid);
            }
          }
        }

        /* minor extents */
        real d1m = color().dys(j);
        real d1c = color().dyc(j);
        real d1p = color().dyn(j);
        real d2m = color().dzb(k);
        real d2c = color().dzc(k);
        real d2p = color().dzt(k);

        /* correct near walls */
        if(wall_indicator[nof-1][nof]<0.5) {
          d1m = d1c;          
        }
        if(wall_indicator[nof+1][nof]<0.5) {
          d1p = d1c;
        }
        if(wall_indicator[nof][nof-1]<0.5) {
          d2m = d2c;
        }
        if(wall_indicator[nof][nof+1]<0.5) {
          d2p = d2c;
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,wall_indicator,
                       mMax,imin,imax,
                       d1m,d1c,d1p,d2m,d2c,d2p,
                       cang_m, cang_p, 
                       max_n,nny,nnz,
                       bflag_struct.jfull,bflag_struct.kfull,
                       color().xc(i),
                       kappa[i][j][k],tempflag[i][j][k]
                       ,i,j,k
                      );

      } else if(mMax==Comp::j()) {
        if(!bflag_struct.jfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int jmin=-mof;
        int jmax= mof;  /* normal stencil size 2*mof+1 */

        /* limit stencil size for cut-stencil */
        if(bflag_struct.jminc) jmin=std::max(-mof,sj()-j);
        if(bflag_struct.jmaxc) jmax=std::min( mof,ej()-j);

        /* because of ghost grid in solid/walls */
        int jmin_grid = jmin;
        int jmax_grid = jmax;
        
        /* update at walls */
        if( bflag_struct.jminw && -mof<sj()-j ) {
          jmin = sj()-j-1;
          jmin_grid = jmin+1;
        }
        if( bflag_struct.jmaxw &&  mof>ej()-j ) {
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

        /* cangle values */
        real cang_m(-10.), cang_p(-10.);

        if(jmax_grid!=jmax) {
          cang_p = cangle(i,j+jmax_grid,k);
        }
        if(jmin_grid!=jmin) {
          cang_m = cangle(i,j+jmin_grid,k);
        }

        /* fill wall_indicator */
        for(int kk(-nof); kk<=nof; ++kk) {
          for(int ii(-nof); ii<=nof; ++ii) {
            /* central position has the contact angle */
            if(kk==0&&ii==0) {
              wall_indicator[kk+nof][ii+nof] = cangle(i,j,k); 
            /* indicator is reset */
            } else if(dom->ibody().off(i+ii,j,k+kk)) {
              wall_indicator[kk+nof][ii+nof] = 0.0;
            } else {
              wall_indicator[kk+nof][ii+nof] = 1.0;
            }
          }
        }

        /* fill stencil */
        for(int kk(-nof); kk<=nof; ++kk) {
          for(int ii(-nof); ii<=nof; ++ii) {
            /* stencil is reset with negative values representing undefined vals */
            for(int jj(-mof); jj<=mof; ++jj) {
              gridstencil[kk+nof][ii+nof][jj+mof] = 0.0;
              stencil[kk+nof][ii+nof][jj+mof] = -1.0;
            }
            for(int jj(jmin); jj<=jmax; ++jj) {
              gridstencil[kk+nof][ii+nof][jj+mof] = color().dyc(j+jj);
              stencil[kk+nof][ii+nof][jj+mof] = std::min(1.0,std::max(0.0,color()[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(jmin_grid != jmin) {
          for(int kk(-nof); kk<=nof; ++kk) {
            for(int ii(-nof); ii<=nof; ++ii) {
              gridstencil[kk+nof][ii+nof][jmin+mof] = color().dyc(j+jmin_grid);
            }
          }
        }
        if(jmax_grid != jmax) {
          for(int kk(-nof); kk<=nof; ++kk) {
            for(int ii(-nof); ii<=nof; ++ii) {
              gridstencil[kk+nof][ii+nof][jmax+mof] = color().dyc(j+jmax_grid);
            }
          }
        }

        /* minor extents */
        real d1m = color().dzb(k);
        real d1c = color().dzc(k);
        real d1p = color().dzt(k);
        real d2m = color().dxw(i);
        real d2c = color().dxc(i);
        real d2p = color().dxe(i);

        /* correct near walls */
        if(wall_indicator[nof-1][nof]<0.5) {
          d1m = d1c;
        }
        if(wall_indicator[nof+1][nof]<0.5) {
          d1p = d1c;
        }
        if(wall_indicator[nof][nof-1]<0.5) {
          d2m = d2c;
        }
        if(wall_indicator[nof][nof+1]<0.5) {
          d2p = d2c;
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,wall_indicator,
                       mMax,jmin,jmax,
                       d1m,d1c,d1p,d2m,d2c,d2p,
                       cang_m, cang_p, 
                       max_n,nnz,nnx,
                       bflag_struct.kfull,bflag_struct.ifull,
                       color().xc(i),
                       kappa[i][j][k],tempflag[i][j][k]
                       ,i,j,k
                      );

      } else { 
        if(!bflag_struct.kfull) {
          boil::oout<<"Curv_HF: Pseudo direction selected at "<<i<<" "<<j<<" "
                    <<k<<" with mMax="<<mMax<<" ; exiting."<<boil::endl;
          exit(0);
        }

        /* check stencil size */
        int kmin=-mof;
        int kmax= mof;  /* normal stencil size 2*mof+1 */ 

        /* limit stencil size for cut-stencil */
        if(bflag_struct.kminc) kmin=std::max(-mof,sk()-k);
        if(bflag_struct.kmaxc) kmax=std::min( mof,ek()-k);

        /* because of ghost grid in solid/walls */
        int kmin_grid = kmin;
        int kmax_grid = kmax;
        
        /* update at walls */
        if( bflag_struct.kminw && -mof<sk()-k ) {
          kmin = sk()-k-1;
          kmin_grid = kmin+1;
        }
        if( bflag_struct.kmaxw &&  mof>ek()-k ) {
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

        /* cangle values */
        real cang_m(-10.), cang_p(-10.);

        if(kmax_grid!=kmax) {
          cang_p = cangle(i,j,k+kmax_grid);
        }
        if(kmin_grid!=kmin) {
          cang_m = cangle(i,j,k+kmin_grid);
        }

        /* fill wall_indicator */
        for(int ii(-nof); ii<=nof; ++ii) {
          for(int jj(-nof); jj<=nof; ++jj) {
            /* central position has the contact angle */
            if(ii==0&&jj==0) {
              wall_indicator[ii+nof][jj+nof] = cangle(i,j,k); 
            /* indicator is reset */
            } else if(dom->ibody().off(i+ii,j+jj,k)) {
              wall_indicator[ii+nof][jj+nof] = 0.0;
            } else {
              wall_indicator[ii+nof][jj+nof] = 1.0;
            }
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
              gridstencil[ii+nof][jj+nof][kk+mof] = color().dzc(k+kk);
              stencil[ii+nof][jj+nof][kk+mof] = std::min(1.0,std::max(0.0,color()[i+ii][j+jj][k+kk]));
            }
          }
        }

        /* correct grid stencil near solid/walls */
        if(kmin_grid != kmin) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int jj(-nof); jj<=nof; ++jj) {
              gridstencil[ii+nof][jj+nof][kmin+mof] = color().dzc(k+kmin_grid);
            }
          }
        }
        if(kmax_grid != kmax) {
          for(int ii(-nof); ii<=nof; ++ii) {
            for(int jj(-nof); jj<=nof; ++jj) {
              gridstencil[ii+nof][jj+nof][kmax+mof] = color().dzc(k+kmax_grid);
            }
          }
        }

        /* minor extents */
        real d1m = color().dxw(i);
        real d1c = color().dxc(i);
        real d1p = color().dxe(i);
        real d2m = color().dys(j);
        real d2c = color().dyc(j);
        real d2p = color().dyn(j);

        /* correct near walls */
        if(wall_indicator[nof-1][nof]<0.5) {
          d1m = d1c; 
        }
        if(wall_indicator[nof+1][nof]<0.5) {
          d1p = d1c;
        }
        if(wall_indicator[nof][nof-1]<0.5) {
          d2m = d2c;
        }
        if(wall_indicator[nof][nof+1]<0.5) {
          d2p = d2c;
        }

        /* calculate curvature */
        curv_HF_kernel(stencil,gridstencil,wall_indicator,
                       mMax,kmin,kmax,
                       d1m,d1c,d1p,d2m,d2c,d2p,
                       cang_m, cang_p, 
                       max_n,nnx,nny,
                       bflag_struct.ifull,bflag_struct.jfull,
                       color().xc(i),
                       kappa[i][j][k],tempflag[i][j][k]
                       ,i,j,k
                      );

      }

    } /* tempflag = 1 */
  } /* for ijk */

  kappa.bnd_update_symmetry(); /* copy on symmetry plane */
  tempflag.bnd_update_symmetry(); /* copy on symmetry plane */
  kappa.exchange();
  tempflag.exchange();

#if 0
  /* near-wall calculations */
  /* if(wall_curv_method==CurvMethod::HFmixedXZ()) {
    insert_bc_curv_HFmixed(color(),Comp::i(),Comp::k(),Sign::neg());
  } else */ if(wall_curv_method==CurvMethod::HFparallelXZ()) {
    insert_bc_curv_HFparallel(color(),Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::HFnormalXZ()) {
    insert_bc_curv_HFnormal(color(),Comp::i(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::HFnormalZ()) {
    insert_bc_curv_HFnormal(color(),Comp::k(),Sign::neg());
  } else if(wall_curv_method==CurvMethod::none()) {
  } else {
    /* default */
    boil::oout<<"VOF::curvHF: Wall curvature calculation method not set properly!"
              <<" Exiting."<<boil::endl;
    exit(0);
  }
#endif

#if 1
  /* enforce curvature in wall adjacent cells = -div.(norm) */
  /* this is a temporal solution because HF doesn't work well in 3D */
  if (wall_curv_method==CurvMethod::DivNorm()) {
    insert_bc_curv_divnorm();
  }
#endif

  /* detect too large kappa */
  real kappa_sphere = 2.0/dom->dxyz_min();
  for_ijk(i,j,k) {
    if((tempflag[i][j][k]==1) && (kappa[i][j][k]>10*kappa_sphere)) {
      boil::aout<<"vof_curv_HF:WARNING!!!  Large curvature: kappa= "<<kappa[i][j][k]
                <<" dxyz_min= "<<dom->dxyz_min()<<" at i= "<<i<<" j= "<<j
                <<" k= "<<k<<" proc= "<<boil::cart.iam()<<"\n";
    }
  }

  /* extrapolate kappa */
  stmp  = kappa;
  tempflag2 = tempflag;
  
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
  /* visualize tempflag */
  if(time->current_step()==1) {
     boil::plot->plot(color(),nx,ny,nz, "clr-nx-ny-nz", time->current_step());
    for_ijk(i,j,k){
      stmp[i][j][k]=tempflag[i][j][k];
    }
    boil::plot->plot(color(),kappa,stmp, "clr-kappa-tempflag", time->current_step());
    exit(0);
  }
#endif

  return;
}
