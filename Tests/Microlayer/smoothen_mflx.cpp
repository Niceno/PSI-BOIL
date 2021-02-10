  /**** smoothen mass flux, assuming 2D ****/
  stmp  = mflx;
  
  /* flagging */
  tempflag = 0;
  for_vijk(tempflag,i,j,k) {
    if(abs((*(conc.topo->iflag))[i][j][k])==1) {
      tempflag[i][j][k] = 1;
    }
  }
  tempflag.bnd_update_symmetry(); /* copy on symmetry plane */
  tempflag.exchange();

  tempflag2 = tempflag;

  /* extend the mflx field */
  for(int iloop=1; iloop<6; iloop++) { 
    for_vijk(tempflag,i,j,k) {
      //if(d.ibody().off(i,j,k)) continue;
      if(tempflag[i][j][k]==0) {
        int inb = std::min(1,tempflag[i-1][j][k]) + std::min(1,tempflag[i+1][j][k])
                + std::min(1,tempflag[i][j][k-1]) + std::min(1,tempflag[i][j][k+1]);
        if(inb >= 1) {
          stmp[i][j][k] = (real(std::min(1,tempflag[i-1][j][k])) * mflx[i-1][j][k]
                        +  real(std::min(1,tempflag[i+1][j][k])) * mflx[i+1][j][k]
                        +  real(std::min(1,tempflag[i][j][k-1])) * mflx[i][j][k-1]
                        +  real(std::min(1,tempflag[i][j][k+1])) * mflx[i][j][k+1])
                        /real(inb);
          tempflag2[i][j][k] = 2;  /* tempflag=2 for extrapolated */
        }
      }
    }
    stmp.bnd_update_symmetry(); /* copy on symmetry plane */
    tempflag2.bnd_update_symmetry(); /* copy on symmetry plane */
    stmp.exchange();
    tempflag2.exchange();
    mflx = stmp;
    tempflag = tempflag2;
  }

  mflx.exchange_all();

  /* smoothen, assuming 2D */
  for(int iloop=1; iloop<2; iloop++) {
    for_vijk(tempflag,i,j,k) {
      //if(d.ibody().off(i,j,k)) continue;
      if(tempflag[i][j][k]>0) {
#if 0 /* Gaussian blur, 3x3 */
        stmp[i][j][k] = (mflx[i  ][j][k]*4.
                      +  mflx[i-1][j][k]*2.
                      +  mflx[i+1][j][k]*2.
                      +  mflx[i][j][k-1]*2.
                      +  mflx[i][j][k+1]*2.
                      +  mflx[i-1][j][k-1]*1.
                      +  mflx[i-1][j][k+1]*1.
                      +  mflx[i+1][j][k-1]*1.
                      +  mflx[i+1][j][k+1]*1.
                      )/16.;
#else /* Gaussian blur, 5x5 */
        stmp[i][j][k] = (mflx[i  ][j][k]*41.
                      +  mflx[i-1][j][k]*26.
                      +  mflx[i+1][j][k]*26.
                      +  mflx[i][j][k-1]*26.
                      +  mflx[i][j][k+1]*26.
                      +  mflx[i-1][j][k-1]*16.
                      +  mflx[i-1][j][k+1]*16.
                      +  mflx[i+1][j][k-1]*16.
                      +  mflx[i+1][j][k+1]*16.
                      +  mflx[i-2][j][k]*7.
                      +  mflx[i+2][j][k]*7.
                      +  mflx[i][j][k-2]*7.
                      +  mflx[i][j][k+2]*7.
                      +  mflx[i-2][j][k-1]*4.
                      +  mflx[i-2][j][k+1]*4.
                      +  mflx[i+2][j][k-1]*4.
                      +  mflx[i+2][j][k+1]*4.
                      +  mflx[i-1][j][k-2]*4.
                      +  mflx[i-1][j][k+2]*4.
                      +  mflx[i+1][j][k-2]*4.
                      +  mflx[i+1][j][k+2]*4.
                      +  mflx[i-2][j][k-2]*1.
                      +  mflx[i-2][j][k+2]*1.
                      +  mflx[i+2][j][k-2]*1.
                      +  mflx[i+2][j][k+2]*1.
                      )/273.;
#endif
      }
    }
    stmp.bnd_update_symmetry(); /* copy on symmetry plane */
    stmp.exchange_all();
    mflx = stmp;
  }
