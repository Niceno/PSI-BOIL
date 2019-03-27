#include "vof.h"

void VOF::superpose() {

  Comp m;

#if 1
  /*----------------------- x-y-z & x-z-y ---------------------*/
  m = Comp::u();
  //boil::oout<<"presup 36 47 48 | "<<stmp[36][47][48]/dV(36,47,48)<<" "<<sosflux[m][36][47][48]/dV(36,47,48)<<" "<<fluxmax[m][36][47][48]/dV(36,47,48)<<" | "<<sosflux[m][37][47][48]/dV(36,47,48)<<" "<<fluxmax[m][37][47][48]/dV(36,47,48)<<boil::endl;

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
    stmp4[iup][j][k] -= fsgn*fval;
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
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in y-direction */
    stmp5[i][jup][k] -= fsgn*fval;
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
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in z-direction */
    stmp6[i][j][kup] -= fsgn*fval;
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
    stmp4[i][jup][k] -= fsgn*fval;
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
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in x-direction */
    stmp5[iup][j][k] -= fsgn*fval;
  }
  stmp5.exchange();

  //boil::oout<<"presup 36 47 48 | "<<stmp[36][47][48]/dV(36,47,48)<<" "<<sosflux[m][36][47][48]/dV(36,47,48)<<" "<<fluxmax[m][36][47][48]/dV(36,47,48)<<" | "<<sosflux[m][37][47][48]/dV(36,47,48)<<" "<<fluxmax[m][37][47][48]/dV(36,47,48)<<boil::endl;
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
    /* update z-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in z-direction */
    stmp6[i][j][kup] -= fsgn*fval;
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
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }
  //boil::oout<<"presup 36 47 48 | "<<stmp[36][47][48]/dV(36,47,48)<<" "<<sosflux[m][36][47][48]/dV(36,47,48)<<" "<<fluxmax[m][36][47][48]/dV(36,47,48)<<" | "<<sosflux[m][37][47][48]/dV(36,47,48)<<" "<<fluxmax[m][37][47][48]/dV(36,47,48)<<boil::endl;

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
    stmp4[i][j][kup] -= fsgn*fval;
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
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in x-direction */
    stmp5[iup][j][k] -= fsgn*fval;
  }
  stmp5.exchange();

  //boil::oout<<"presup 36 47 48 | "<<stmp[36][47][48]/dV(36,47,48)<<" "<<sosflux[m][36][47][48]/dV(36,47,48)<<" "<<fluxmax[m][36][47][48]/dV(36,47,48)<<" | "<<sosflux[m][37][47][48]/dV(36,47,48)<<" "<<fluxmax[m][37][47][48]/dV(36,47,48)<<boil::endl;
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
    /* update y-flux */
    sosflux[m][i][j][k] += fval/6.0;
    /* update stmp in y-direction */
    stmp6[i][jup][k] -= fsgn*fval;
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
    /* update x-flux */
    sosflux[m][i][j][k] += fval/6.0;
  }
  //boil::oout<<"presup 36 47 48 | "<<stmp[36][47][48]/dV(36,47,48)<<" "<<sosflux[m][36][47][48]/dV(36,47,48)<<" "<<fluxmax[m][36][47][48]/dV(36,47,48)<<" | "<<sosflux[m][37][47][48]/dV(36,47,48)<<" "<<fluxmax[m][37][47][48]/dV(36,47,48)<<boil::endl;
#endif

  /* update stmp */
  sosflux.exchange();

#if 0 
    if(time->current_step()==2 ){
    m = Comp::u();
    for_wvmi(sosflux,m,i)
       boil::oout<<i<<" "<<stmp[i][48][48]<<" "<<sosflux[m][i][48][48]<<" "<<fluxmax[m][i][48][48]
                 <<" | "<<stmp[i][48][49]<<" "<<sosflux[m][i][48][49]<<" "<<fluxmax[m][i][48][49]
                 <<" | "<<stmp[i][49][48]<<" "<<sosflux[m][i][49][48]<<" "<<fluxmax[m][i][49][48]
                 <<" | "<<stmp[i][49][49]<<" "<<sosflux[m][i][49][49]<<" "<<fluxmax[m][i][49][49]<<boil::endl;
    boil::plot->plot(sosflux,stmp,"flux-color",time->current_step());
    exit(0);
    }
#endif

#if 1 
  m = Comp::u();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i-1][j][k] = stmp[i-1][j][k] - sosflux[m][i][j][k];
    stmp[i  ][j][k] = stmp[i  ][j][k] + sosflux[m][i][j][k];

    //if(stmp[i-1][j][k]<.0) boil::oout<<"xm: "<<i-1<<" "<<j<<" "<<k<<" | "<<stmp[i-1][j][k]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
    //if(stmp[i  ][j][k]<.0) boil::oout<<"xp: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i  ][j][k]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
  }  
  m = Comp::v();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j-1][k] = stmp[i][j-1][k] - sosflux[m][i][j][k];
    stmp[i][j  ][k] = stmp[i][j  ][k] + sosflux[m][i][j][k];

    //if(stmp[i][j-1][k]<.0) boil::oout<<"ym: "<<i<<" "<<j-1<<" "<<k<<" | "<<stmp[i][j-1][k]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
    //if(stmp[i][j  ][k]<.0) boil::oout<<"yp: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j  ][k]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
  }
  m = Comp::w();
  for_wvmijk(sosflux,m,i,j,k) {
    stmp[i][j][k-1] = stmp[i][j][k-1] - sosflux[m][i][j][k];
    stmp[i][j][k  ] = stmp[i][j][k  ] + sosflux[m][i][j][k];

    //if(stmp[i][j][k-1]<.0) boil::oout<<"zm: "<<i<<" "<<j<<" "<<k-1<<" | "<<stmp[i][j][k-1]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
    //if(stmp[i][j][k  ]<.0) boil::oout<<"zp: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j][k  ]/dV(i,j,k)<<" "<<sosflux[m][i][j][k]/dV(i,j,k)<<" "<<fluxmax[m][i][j][k]/dV(i,j,k)<<boil::endl;
  }  

  for_ijk(i,j,k) {
    if(stmp[i][j][k]/dV(i,j,k)>1.0+boil::micro) {
      boil::oout<<"superpose: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j][k]/dV(i,j,k)<<" | "<<(*u)[Comp::u()][i][j][k]<<" "<<(*u)[Comp::u()][i+1][j][k]<<" | "<<(*u)[Comp::v()][i][j][k]<<" "<<(*u)[Comp::v()][i][j+1][k]<<" | "<<(*u)[Comp::w()][i][j][k]<<" "<<(*u)[Comp::w()][i][j][k+1]<<" | "<<sosflux[Comp::u()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::u()][i+1][j][k]/dV(i,j,k)<<" | "<<sosflux[Comp::v()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::v()][i][j+1][k]/dV(i,j,k)<<" | "<<sosflux[Comp::w()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::w()][i][j][k+1]/dV(i,j,k)<<" | "<<uliq[Comp::u()][i][j][k]<<" "<<uliq[Comp::u()][i+1][j][k]<<" | "<<uliq[Comp::v()][i][j][k]<<" "<<uliq[Comp::v()][i][j+1][k]<<" | "<<uliq[Comp::w()][i][j][k]<<" "<<uliq[Comp::w()][i][j][k+1]<<boil::endl;
     boil::oout<<boil::endl;
   }
 }
#else
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

  for_ijk(i,j,k) {
    if(stmp[i][j][k]/dV(i,j,k)>1.0+boil::micro) {
      boil::oout<<"superpose: "<<i<<" "<<j<<" "<<k<<" | "<<stmp[i][j][k]/dV(i,j,k)<<" | "<<(*u)[Comp::u()][i][j][k]<<" "<<(*u)[Comp::u()][i+1][j][k]<<" | "<<(*u)[Comp::v()][i][j][k]<<" "<<(*u)[Comp::v()][i][j+1][k]<<" | "<<(*u)[Comp::w()][i][j][k]<<" "<<(*u)[Comp::w()][i][j][k+1]<<" | "<<sosflux[Comp::u()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::u()][i+1][j][k]/dV(i,j,k)<<" | "<<sosflux[Comp::v()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::v()][i][j+1][k]/dV(i,j,k)<<" | "<<sosflux[Comp::w()][i][j][k]/dV(i,j,k)<<" "<<sosflux[Comp::w()][i][j][k+1]/dV(i,j,k)<<" | "<<uliq[Comp::u()][i][j][k]<<" "<<uliq[Comp::u()][i+1][j][k]<<" | "<<uliq[Comp::v()][i][j][k]<<" "<<uliq[Comp::v()][i][j+1][k]<<" | "<<uliq[Comp::w()][i][j][k]<<" "<<uliq[Comp::w()][i][j][k+1]<<boil::endl;
     boil::oout<<boil::endl;
   }
  } 
#endif

  return;
}
