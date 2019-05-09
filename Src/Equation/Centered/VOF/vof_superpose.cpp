#include "vof.h"
//#define SYMMETRIC
//#define HYPER

void VOF::superpose() {

#if 0
  /* reset */
  for_m(m) 
    for_avmijk(sosflux,m,i,j,k) {
      sosflux[m][i][j][k] = 0.0;
    }

  Comp m;
#if 1
  /*----------------------- x-y-z & x-z-y ---------------------*/
  m = Comp::u();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp4[i][j][k] = stmp[i][j][k];
  }

  /* update stmp in x-direction */
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(fval<0.0) {
      iup = i; 
      idn = i-1;
    }
#ifdef HYPER
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[iup][j][k])); 
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(idn,j,k)-stmp4[idn][j][k])); 
#endif
    /* update x-flux */
    sosflux[m][i][j][k] += fval/3.0;
    /* update stmp in x-direction */
    stmp4[iup][j][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp4[idn][j][k] += fsgn*fval;
#endif
  }  
  stmp4.exchange();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp5[i][j][k] = stmp4[i][j][k];
    stmp6[i][j][k] = stmp4[i][j][k];
  }

  /* superpose in y-direction */
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(fval<0.0) {
      jup = j; 
      jdn = j-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[i][jup][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,jdn,k)-stmp4[i][jdn][k])); 
#endif
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in y-direction */
    stmp5[i][jup][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp5[i][jdn][k] += fsgn*fval;
#endif
  }
  stmp5.exchange();

  /* superpose in z-direction */
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(fval<0.0) {
      kup = k; 
      kdn = k-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(stmp4[i][j][kup],0.0)); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,j,kdn)-stmp4[i][j][kdn])); 
#endif
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in z-direction */
    stmp6[i][j][kup] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp6[i][j][kdn] += fsgn*fval;
#endif
  }
  stmp6.exchange();

  /* superpose in yz-direction */
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(fval<0.0) {
      kup = k; 
      kdn = k-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp5[i][j][kup])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,j,kdn)-stmp5[i][j][kdn])); 
#endif
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }

  /* superpose in zy-direction */
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(fval<0.0) {
      jup = j; 
      jdn = j-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp6[i][jup][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,jdn,k)-stmp6[i][jdn][k])); 
#endif
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }

  /*----------------------- y-x-z & y-z-x ---------------------*/
  m = Comp::v();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp4[i][j][k] = stmp[i][j][k];
  }

  /* update stmp in y-direction */
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(fval<0.0) {
      jup = j; 
      jdn = j-1;
    }
#ifdef HYPER
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[i][jup][k])); 
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,jdn,k)-stmp4[i][jdn][k])); 
#endif
    /* update y-flux */
    sosflux[m][i][j][k] += fval/3.0;
    /* update stmp in y-direction */
    stmp4[i][jup][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp4[i][jdn][k] += fsgn*fval;
#endif
  }  
  stmp4.exchange();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp5[i][j][k] = stmp4[i][j][k];
    stmp6[i][j][k] = stmp4[i][j][k];
  }

  /* superpose in x-direction */
  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(fval<0.0) {
      iup = i; 
      idn = i-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[iup][j][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(idn,j,k)-stmp4[idn][j][k])); 
#endif
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in x-direction */
    stmp5[iup][j][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp5[idn][j][k] += fsgn*fval;
#endif
  }
  stmp5.exchange();

  /* superpose in z-direction */
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(fval<0.0) {
      kup = k; 
      kdn = k-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[i][j][kup])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,j,kdn)-stmp4[i][j][kdn])); 
#endif
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in z-direction */
    stmp6[i][j][kup] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp6[i][j][kdn] += fsgn*fval;
#endif
  }
  stmp6.exchange();

  /* superpose in xz-direction */
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(fval<0.0) {
      kup = k; 
      kdn = k-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp5[i][j][kup])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,j,kdn)-stmp5[i][j][kdn])); 
#endif
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }

  /* superpose in zx-direction */
  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(fval<0.0) {
      iup = i; 
      idn = i-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp6[iup][j][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(idn,j,k)-stmp6[idn][j][k])); 
#endif
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }

  /*----------------------- z-x-y & z-y-x ---------------------*/
  m = Comp::w();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp4[i][j][k] = stmp[i][j][k];
  }

  /* update stmp in z-direction */
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind k-index */
    int kup(k-1), kdn(k);
    if(fval<0.0) {
      kup = k; 
      kdn = k-1;
    }
#ifdef HYPER
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[i][j][kup])); 
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,j,kdn)-stmp4[i][j][kdn])); 
#endif
    /* update z-flux */
    sosflux[m][i][j][k] += fval/3.0;
    /* update stmp in z-direction */
    stmp4[i][j][kup] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp4[i][j][kdn] += fsgn*fval;
#endif
  }  
  stmp4.exchange();

  /* reset */
  for_avijk(stmp4,i,j,k) {
    stmp5[i][j][k] = stmp4[i][j][k];
    stmp6[i][j][k] = stmp4[i][j][k];
  }

  /* superpose in x-direction */
  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(fval<0.0) {
      iup = i; 
      idn = i-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[iup][j][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(idn,j,k)-stmp4[idn][j][k])); 
#endif
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in x-direction */
    stmp5[iup][j][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp5[idn][j][k] += fsgn*fval;
#endif
  }
  stmp5.exchange();

  /* superpose in y-direction */
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(fval<0.0) {
      jup = j; 
      jdn = j-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp4[i][jup][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,jdn,k)-stmp4[i][jdn][k])); 
#endif
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in y-direction */
    stmp6[i][jup][k] -= fsgn*fval;
#ifdef SYMMETRIC
    stmp6[i][jdn][k] += fsgn*fval;
#endif
  }
  stmp6.exchange();

  /* superpose in xy-direction */
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind j-index */
    int jup(j-1), jdn(j);
    if(fval<0.0) {
      jup = j; 
      jdn = j-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp5[i][jup][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(i,jdn,k)-stmp5[i][jdn][k])); 
#endif
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }

  /* superpose in yx-direction */
  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    real fval = fluxmax[m][i][j][k];
    real fsgn = (fval>0.0) - (fval<0.0);
    /* upwind i-index */
    int iup(i-1), idn(i);
    if(fval<0.0) {
      iup = i; 
      idn = i-1;
    }
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,stmp6[iup][j][k])); 
#ifdef SYMMETRIC
    fval = fsgn*std::min(fsgn*fval,std::max(0.0,dV(idn,j,k)-stmp6[idn][j][k])); 
#endif
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }
#endif

  /* update stmp */
  sosflux.exchange();

  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i-1][j][k] = stmp[i-1][j][k] - sosflux[m][i][j][k];
    stmp[i  ][j][k] = stmp[i  ][j][k] + sosflux[m][i][j][k];
  }  
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j-1][k] = stmp[i][j-1][k] - sosflux[m][i][j][k];
    stmp[i][j  ][k] = stmp[i][j  ][k] + sosflux[m][i][j][k];
  }
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j][k-1] = stmp[i][j][k-1] - sosflux[m][i][j][k];
    stmp[i][j][k  ] = stmp[i][j][k  ] + sosflux[m][i][j][k];
  }  

  #if 0
  for_ijk(i,j,k) {
    if(stmp[i][j][k]/dV(i,j,k)>1.0+boil::micro) {
      boil::oout<<"superpose: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j][k]/dV(i,j,k)<<" | "<<(*u)[Comp::u()][i][j][k]<<" "<<(*u)[Comp::u()][i+1][j][k]<<" | "<<(*u)[Comp::v()][i][j][k]<<" "<<(*u)[Comp::v()][i][j+1][k]<<" | "<<(*u)[Comp::w()][i][j][k]<<" "<<(*u)[Comp::w()][i][j][k+1]<<" | "<<sosflux[Comp::u()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::u()][i+1][j][k]/dV(i,j,k)<<" | "<<sosflux[Comp::v()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::v()][i][j+1][k]/dV(i,j,k)<<" | "<<sosflux[Comp::w()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::w()][i][j][k+1]/dV(i,j,k)<<" | "<<uliq[Comp::u()][i][j][k]<<" "<<uliq[Comp::u()][i+1][j][k]<<" | "<<uliq[Comp::v()][i][j][k]<<" "<<uliq[Comp::v()][i][j+1][k]<<" | "<<uliq[Comp::w()][i][j][k]<<" "<<uliq[Comp::w()][i][j][k+1]<<boil::endl;
     boil::oout<<boil::endl;
   }
 }
  #endif
#else
  Comp m;

  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i-1][j][k] = stmp[i-1][j][k] - fluxmax[m][i][j][k];
    stmp[i  ][j][k] = stmp[i  ][j][k] + fluxmax[m][i][j][k];
  }  
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j-1][k] = stmp[i][j-1][k] - fluxmax[m][i][j][k];
    stmp[i][j  ][k] = stmp[i][j  ][k] + fluxmax[m][i][j][k];
  }
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j][k-1] = stmp[i][j][k-1] - fluxmax[m][i][j][k];
    stmp[i][j][k  ] = stmp[i][j][k  ] + fluxmax[m][i][j][k];
  }  

  #if 0
  for_ijk(i,j,k) {
    if(stmp[i][j][k]/dV(i,j,k)>1.0+boil::micro) {
      boil::oout<<"superpose: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j][k]/dV(i,j,k)<<" | "<<(*u)[Comp::u()][i][j][k]<<" "<<(*u)[Comp::u()][i+1][j][k]<<" | "<<(*u)[Comp::v()][i][j][k]<<" "<<(*u)[Comp::v()][i][j+1][k]<<" | "<<(*u)[Comp::w()][i][j][k]<<" "<<(*u)[Comp::w()][i][j][k+1]<<" | "<<sosflux[Comp::u()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::u()][i+1][j][k]/dV(i,j,k)<<" | "<<sosflux[Comp::v()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::v()][i][j+1][k]/dV(i,j,k)<<" | "<<sosflux[Comp::w()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::w()][i][j][k+1]/dV(i,j,k)<<" | "<<uliq[Comp::u()][i][j][k]<<" "<<uliq[Comp::u()][i+1][j][k]<<" | "<<uliq[Comp::v()][i][j][k]<<" "<<uliq[Comp::v()][i][j+1][k]<<" | "<<uliq[Comp::w()][i][j][k]<<" "<<uliq[Comp::w()][i][j][k+1]<<boil::endl;
     boil::oout<<boil::endl;
   }
  } 
  #endif
#endif

  return;
}
