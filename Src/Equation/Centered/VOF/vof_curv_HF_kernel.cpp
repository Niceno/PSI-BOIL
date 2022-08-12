#include "vof.h"
#define LOPEZ_COLORING
#define LOPEZ_SMOOTHING

/*-------------------------------+
|  curvature calculation kernel  |
+-------------------------------*/
/* inputs:
   - stencil of color values
   - stencil of grid spacings in the dominant direction
   - stencil of wall indicators
   - direction of construction
   - starting and ending indices in the dominant direction
   - grid spacings in the first minorext direction
   - grid spacings in the second minorext direction
   - contact angles at bottom and top
   - normal vector component in the dominant direction, pointing to gas
   - bool indicator if first minorext direction is a non-pseudo direction
   - bool indicator if second minorext direction is a non-pseudo direction
   - radial distance from origin (used in axisym)
   
   outputs:
   - kappa value (if computed)
   - tempflag value (if kappa not computed)

   debug:
   - i,j,k coordinates
*/
void VOF::curv_HF_kernel(arr3D & stencil, const arr3D & gridstencil,
                         const arr2D & wall_indicator, const Comp & mcomp,
                         const int imin, const int imax,
                         const real d1m, const real d1c, const real d1p,
                         const real d2m, const real d2c, const real d2p,
                         const real cng_m, const real cng_p,
                         const real max_n, const real nn_j, const real nn_k,
                         const bool truedir1, const bool truedir2,
                         const real xcent,
                         real & kap, int & flag
                         ,const int i, const int j, const int k
                         ) const {

#if 0
  if (k==4+sk()&&j==2+sj()) {
    std::cout<<"HF_kernel: "<<i<<" "<<j<<" "<<k<<" "<<"\n";
    std::cout<<"cng_m "<<cng_m<<" cng_p "<<cng_p<<" max_n= "<<max_n<<" nn_j "<<nn_j<<" nn_k "<<nn_k<<" imin "<<imin<<" imax "<<imax<<"\n";
  }
#endif
#if 0
  if (k==4+sk()&&j==2+sj()&&i==117) {
    std::cout<<"curv_HF_kernel:stencil: "<<stencil[1][1][2]<<" "<<stencil[1][1][3]<<" "
    <<stencil[1][1][4]<<" "<<stencil[1][1][5]<<" "<<stencil[1][1][6]<<"\n";
    std::cout<<"truedir1 "<<truedir1<<" truedir2 "<<truedir2<<"\n";
  }
#endif

  const int & mof = hf_set.mof;
  const int & nof = hf_set.nof;

  real cang_m = cng_m;
  real cang_p = cng_p;

  real nnorm = nn_j*nn_j+nn_k*nn_k;
  real nnj = nn_j/sqrt(nnorm+boil::pico);
  real nnk = nn_k/sqrt(nnorm+boil::pico);

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
    /* also contact angles must be inverted */
    if(cang_m>-1.)
      cang_m = boil::pi-cang_m;
    if(cang_p>-1.)
      cang_p = boil::pi-cang_p;
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

#if 0
  if (k==4+sk()&&j==2+sj()&&i==117) {
    std::cout<<"curv_HF_kernel:stencil-invert: "<<stencil[1][1][2]<<" "<<stencil[1][1][3]<<" "
    <<stencil[1][1][4]<<" "<<stencil[1][1][5]<<" "<<stencil[1][1][6]<<"\n";
    std::cout<<"cang_m "<<cang_m<<" cang_p "<<cang_p<<"\n";
  }
#endif

  /* calculate hc_limit */
  real nhc_0(0.),nhc_1(gridstencil[nof][nof][mof]);
  for(int ii(imin); ii<0; ++ii) {
    nhc_0 += gridstencil[nof][nof][ii+mof];
    nhc_1 += gridstencil[nof][nof][ii+mof];
  }

  /* height stencil */
  arr2D heights;
  heights.resize(hf_set.minorext);
  for(auto & hs : heights) {
    hs.resize(hf_set.minorext);
    for(auto & h : hs) 
      h = 0.;
  }

  /* fill the height stencil */
  for(int jj(-nof); jj<=nof; ++jj) {
    for(int kk(-nof); kk<=nof; ++kk) {
      for(int ii(imin); ii<=imax; ++ii) {
        heights[jj+nof][kk+nof] +=     stencil[jj+nof][kk+nof][ii+mof]
                                  *gridstencil[jj+nof][kk+nof][ii+mof];
      }
    }
  }

#if 0
  if (k==4+sk()&&j==2+sj()&&i==117) {
    std::cout<<"Heights-0: ";
    for(int jj(-nof); jj<=nof; ++jj) {
      for(int kk(-nof); kk<=nof; ++kk) {
        std::cout<<jj+nof<<" "<<kk+nof<<" "<<heights[jj+nof][kk+nof]<<" | ";
      }
    std::cout<<" || ";
    }
    std::cout<<"\n";
  }
#endif
#if 0
  boil::oout<<"heights: ";
  for(int jj(-nof); jj<=nof; ++jj) {
    for(int kk(-nof); kk<=nof; ++kk) {
      boil::oout<<jj+nof<<" "<<kk+nof<<" "<<heights[jj+nof][kk+nof]<<" | ";
    }
    boil::oout<<" || ";
  }
  boil::oout<<boil::endl;
#endif

  /* a slope in any given direction should take into account the
     projection onto the boundary plane, calculating height drop as
                   CABS * < (nx,ny) , (Dx,Dy) >,
     where the former is a contact-angle-based slope (CABS) and the
     latter term represents a scalar product of the normal
     vector projection onto the boundary plane and a distance vector
     in the chosen direction, n cdot D. Note that:
                   n cdot D = ||n|| * ||d|| * cos(alpha),
      where ||n|| = 1 and alpha is the angle between the normal
      vector and the distance vector. Thus, n cdot D is numerically
      equal to the length of the projection of D in the n-dir */

  /* contact-angle correction in major direction */
  if(cang_m>-1.) {
    /* update at walls, otherwise should be zero */
    real height_min = gridstencil[nof][nof][imin+mof];
	
	/* avoid tangens of pi/2 */
	real tancang_m = signum(1.0,sin(cang_m))*boil::unreal;
	real signum_cos_m = signum(1.0,cos(cang_m));
	if(std::abs(signum_cos_m)>0.5) {
	  tancang_m *= signum_cos_m;
	}
	if(std::abs(cang_m-boil::pi/2.)>boil::pico) {
	  tancang_m = tan(cang_m);  
	}	

#if 0
      if (k==4+sk()&&j==2+sj()&&i==117) {
        std::cout<<"tancang_m= "<<tancang_m<<" nof= "<<nof<<"\n";
      }
#endif

    for(int jj(-nof); jj<=nof; ++jj) {
      for(int kk(-nof); kk<=nof; ++kk) {
        if(jj==0&&kk==0)
          continue;
#if 0
        if (k==4+sk()&&j==2+sj()&&i==117) {
          std::cout<<"heights= "<<heights[jj+nof][kk+nof]<<" "<<height_min<<" "<<jj<<" "<<kk<<"\n";
        }
#endif
        if(heights[jj+nof][kk+nof]<height_min) {
          /* only approximate for nof>1 */
          real distj = (jj>0) ? real(jj)*(d1p+d1c)/2. : real(jj)*(d1m+d1c)/2.;
          real distk = (kk>0) ? real(kk)*(d2p+d2c)/2. : real(kk)*(d2m+d2c)/2.;
          
          distj *= nnj*mult;
          distk *= nnk*mult;

          /* we disallow extrapolation in the negative sense, this would
             indicate corrupted normal vector */
          real diff = std::max(0., (distj+distk)*tancang_m );
          heights[jj+nof][kk+nof] = heights[nof][nof]-diff;
        }
      }
    }
  }

  if(cang_p>-1.) {
    /* update at walls, otherwise should be zero */
    real height_max = 0.;
	
	/* avoid tangens of pi/2 */
	real tancang_p = signum(1.0,sin(cang_p))*boil::unreal;
	real signum_cos_p = signum(1.0,cos(cang_p));
	if(std::abs(signum_cos_p)>0.5) {
	  tancang_p *= signum_cos_p;
	}
	if(std::abs(cang_p-boil::pi/2.)>boil::pico) {
	  tancang_p = tan(cang_p);  
	}
	
    for(int ii(imin); ii<imax; ++ii) {
      height_max += gridstencil[nof][nof][ii+mof];
    }
    for(int jj(-nof); jj<=nof; ++jj) {
      for(int kk(-nof); kk<=nof; ++kk) {
        if(jj==0&&kk==0)
          continue;
        if(heights[jj+nof][kk+nof]>height_max) {
          /* approximate for nof>1 */
          real distj = (jj>0) ? real(jj)*(d1p+d1c)/2. : real(jj)*(d1m+d1c)/2.;
          real distk = (kk>0) ? real(kk)*(d2p+d2c)/2. : real(kk)*(d2m+d2c)/2.;

          distj *= nnj*mult;
          distk *= nnk*mult;

          /* we disallow extrapolation in the positive sense, this would
             indicate corrupted normal vector */
          /* the tangens is negative, as well as the dist! */
          real diff = std::max(0., (distj+distk)*tancang_p );
          heights[jj+nof][kk+nof] = heights[nof][nof]+diff;
        }
      }
    }
  }

  /* inactive columns, wall_indicator[nof][nof] is the cangle */
  real cang_c = wall_indicator[nof][nof];
  if(mult<0.) {
    cang_c = boil::pi-cang_c;
  }

  /* avoid cotangens of 0 */
  real cotancang_c = signum(1.0,cos(cang_c))*boil::unreal;
  real signum_sin_c = signum(1.0,sin(cang_c));
  if(std::abs(signum_sin_c)>0.5) {
	cotancang_c *= signum_sin_c;
  }
  if(std::abs(cang_c)>boil::pico && std::abs(cang_c-boil::pi)>boil::pico) {
	cotancang_c = 1./tan(cang_c);  
  }

  for(int jj(-nof); jj<=nof; ++jj) {
    for(int kk(-nof); kk<=nof; ++kk) {
      if(jj==0&&kk==0)
          continue;
      if(wall_indicator[jj+nof][kk+nof]<0.5) {
        /* approximate for nof>1 */
        real distj = (jj>0) ? real(jj)*(d1p+d1c)/2. : real(jj)*(d1m+d1c)/2.;
        real distk = (kk>0) ? real(kk)*(d2p+d2c)/2. : real(kk)*(d2m+d2c)/2.;

        distj *= nnj*mult;
        distk *= nnk*mult;
        real dist = distj+distk;

        heights[jj+nof][kk+nof] = heights[nof][nof] + fabs(dist)*cotancang_c;
      }
    }
  }

#if 0
  if (k==4+sk()&&j==2+sj()&&i==117) {
    std::cout<<"heghts10: ";
    for(int jj(-nof); jj<=nof; ++jj) {
      for(int kk(-nof); kk<=nof; ++kk) {
        std::cout<<jj+nof<<" "<<kk+nof<<" "<<heights[jj+nof][kk+nof]<<" | ";
      }
    std::cout<<" || ";
    }
    std::cout<<"\n";
    exit(0);
  }
#endif
  //if(nhc_limit<nhc && nhc<=(nhc_limit+1.0)) {
  if(nhc_0<heights[nof][nof] && heights[nof][nof]<=nhc_1) {

    kap = calculate_curvature_HF(heights,
                                 d1m, d1c, d1p,
                                 d2m, d2c, d2p,
                                 truedir1, truedir2,
                                 mult, max_n, xcent); 
  } else {
    flag = 0;
  }

#if 0
    if(time->current_step()==305) {
      if(i==25&&j==17&&k==4) {
        std::cout<<"curv_HF_kernel:999:kap= "<<kap<<" "<<flag<<"\n";
      }
    }
#endif

  return;
} 

real VOF::calculate_curvature_HF(const arr2D & heights,
                                 const real d1m, const real d1c, const real d1p,
                                 const real d2m, const real d2c, const real d2p,
                                 const bool truedir1, const bool truedir2,
                                 const real mult, const real max_n,
                                 const real xcent) const {

  const real & hmm = heights[0][0];
  const real & hcm = heights[1][0];
  const real & hpm = heights[2][0];

  const real & hmc = heights[0][1];
  const real & hcc = heights[1][1];
  const real & hpc = heights[2][1];

  const real & hmp = heights[0][2];
  const real & hcp = heights[1][2];
  const real & hpp = heights[2][2];

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
