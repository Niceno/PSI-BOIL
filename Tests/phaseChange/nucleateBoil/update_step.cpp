void update_step(const Scalar & c, Scalar & step, Scalar & sflag){
#if 1
  const real phisurf=0.5;
  for_avijk(sflag,i,j,k) {
    sflag[i][j][k]=0;
  }
  for_vijk(c,i,j,k) {
    if(c[i][j][k]>=phisurf)
      sflag[i][j][k]=1;
  }
  /* i-direction */
  for(int i=c.si()-1; i<=c.ei(); i++){
    for_vjk(c,j,k){
       if((c[i][j][k]-phisurf)*(c[i+1][j][k]-phisurf)<=0.0){
          sflag[i  ][j][k]=2;
          sflag[i+1][j][k]=2;
       }
    }
  }
  /* j-direction */
  for(int j=c.sj()-1; j<=c.ej(); j++){
    for_vik(c,i,k){
      if((c[i][j][k]-phisurf)*(c[i][j+1][k]-phisurf)<=0.0){
          sflag[i][j  ][k]=2;
          sflag[i][j+1][k]=2;
       }
    }
  }
  /* k-direction */
  for(int k=c.sk()-1; k<=c.ek(); k++){
    for_vij(c,i,j){
       if((c[i][j][k]-phisurf)*(c[i][j][k+1]-phisurf)<=0.0){
          sflag[i][j][k  ]=2;
          sflag[i][j][k+1]=2;
       }
    }
  }
  sflag.exchange_all();
  for_avijk(c,i,j,k){
    if(sflag[i][j][k]==2){
      step[i][j][k]=c[i][j][k];
    } else {
      if(c[i][j][k]<0.5){
        step[i][j][k]=0.0;
      } else {
        step[i][j][k]=1.0;
      }
    }
  }
#else
  for_avijk(step,i,j,k){
    step[i][j][k]=c[i][j][k];
  }
#endif
}

