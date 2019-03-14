#include "vof.h"
/******************************************************************************/
#if 0
void VOF::curv_HF_ext(){
  const int ist=-2;
  const int ied=+2;
  int d_i=0, d_j=0, d_k=0;
#if 0
  for_aijk(i,j,k) {
    stmp2[i][j][k]=real(iflag[i][j][k]);
  }
  for_aijk(i,j,k) {
    stmp3[i][j][k]=real(iflagx[i][j][k]);
  }

  if(time->current_step() == 1) {
      boil::plot->plot(phi,stmp2,stmp3, "phi-iflag-iflagx", time->current_step());
  }
#endif
  for_ijk(i,j,k){
    iflag[i][j][k] = iflagx[i][j][k];
  }

  for_ijk(i,j,k){
    if(iflag[i][j][k]==1 && iflagx[i][j][k]!=1){  
      for(int ii=ist; ii<=ied; ii++){
        if(iflagx[i+ii][j][k]==1){  //The algorithm doesn't consider if two interfaces appear in one direction, this can be improved if necessary 
          d_i = ii;
        }  
      }
      for(int jj=ist; jj<=ied; jj++){
        if(iflagx[i][j+jj][k]==1){
          d_j = jj;
        }
      }
      for(int kk=ist; kk<=ied; kk++){
        if(iflagx[i][j][k+kk]==1){
          d_k = kk;
        }
      }
     
      if(d_i==0){
        d_i = 10;
      }
      if(d_j==0){
        d_j = 10;
      } 
      if(d_k==0){
        d_k = 10;
      } 

      /* d_i==d_j==d_k */
      if(fabs(d_i)==fabs(d_j) && fabs(d_j)==fabs(d_k)){
        kappa[i][j][k] = (kappa[i+d_i][j][k]+kappa[i][j+d_j][k]+kappa[i][j][k+d_k])/3.0;
      }   
      /* d_i==d_j>d_k */
      if(fabs(d_i)==fabs(d_j) && fabs(d_j)>fabs(d_k)){
        kappa[i][j][k] = kappa[i][j][k+d_k];
      } 
      /* d_k>d_i==d_j */
      if(fabs(d_k)>fabs(d_i) && fabs(d_i)==fabs(d_j)){
        kappa[i][j][k] = (kappa[i+d_i][j][k]+kappa[i][j+d_j][k])/2.0;
      }
      /* d_i>d_j=d_k */
      if(fabs(d_i)>fabs(d_j) && fabs(d_j)==fabs(d_k)){
        kappa[i][j][k] = (kappa[i][j+d_j][k]+kappa[i][j][k+d_k])/2.0;
      }
      /* d_j==d_k>d_i */
      if(fabs(d_j)==fabs(d_k) && fabs(d_k)>fabs(d_i)){
        kappa[i][j][k] = kappa[i+d_i][j][k];
      }
      /* d_j>d_i=d_k */
      if(fabs(d_j)>fabs(d_i) && fabs(d_i)==fabs(d_k)){
        kappa[i][j][k] = (kappa[i+d_i][j][k]+ kappa[i][j][k+d_k])/2.0;
      }
      /* d_i==d_k>d_j */
      if(fabs(d_i)==fabs(d_k) && fabs(d_k)>fabs(d_j)){
        kappa[i][j][k] = kappa[i][j+d_j][k];
      }
      /* d_i>d_j>d_k */
      if(fabs(d_i)>fabs(d_j) && fabs(d_j)>fabs(d_k)){
        kappa[i][j][k] = kappa[i][j][k+d_k];
      }
      /* d_i>d_k>d_j */
      if(fabs(d_i)>fabs(d_k) && fabs(d_k)>fabs(d_j)){
        kappa[i][j][k] = kappa[i][j+d_j][k];
      }
      /* d_j>d_i>d_k */
      if(fabs(d_j)>fabs(d_i) && fabs(d_i)>fabs(d_k)){
        kappa[i][j][k] = kappa[i][j][k+d_k];
      }
      /* d_j>d_k>d_i */
      if(fabs(d_j)>fabs(d_k) && fabs(d_k)>fabs(d_i)){
        kappa[i][j][k] = kappa[i+d_i][j][k];
      }
      /* d_k>d_i>d_j */
      if(fabs(d_k)>fabs(d_i) && fabs(d_i)>fabs(d_j)){
        kappa[i][j][k] = kappa[i][j+d_j][k];
      }
      /* d_k>d_j>d_i */
      if(fabs(d_k)>fabs(d_j) && fabs(d_j)>fabs(d_i)){
        kappa[i][j][k] = kappa[i+d_i][j][k];
      }
      d_i = 0;
      d_j = 0;
      d_k = 0;
    }
  }
  return;
}
#endif

#if 0
void VOF::curv_HF_ext(){

  int sum=0;
  real flag=0.0;
  real kappa_sum=0.0;

  for_aijk(i,j,k) {
    stmp2[i][j][k]=iflag[i][j][k];
    stmp3[i][j][k]=iflagx[i][j][k];
  }


  for_aijk(i,j,k){
      iflag[i][j][k] = iflagx[i][j][k];
  }
  iflag.exchange();

  for(int layer=1; layer<=3; layer++){
    for_ijk(i,j,k){
      if(iflag[i][j][k]==0){
        for(int ii=-1; ii<=1; ii++){
          sum += iflag[i+ii][j][k];
        } 
        for(int jj=-1; jj<=1; jj++){
          sum += iflag[i][j+jj][k];
        }
        for(int kk=-1; kk<=1; kk++){
          sum += iflag[i][j][k+kk];  
        }
             
        if(sum>0){
          for(int ii=-1; ii<=1; ii++){
            if(iflag[i+ii][j][k]==1){
              kappa_sum += kappa[i+ii][j][k];
              flag += 1.0; 
            }
          } 
          for(int jj=-1; jj<=1; jj++){
            if(iflag[i][j+jj][k]==1){
              kappa_sum += kappa[i][j+jj][k];
              flag += 1.0;
            }
          } 
          for(int kk=-1; kk<=1; kk++){
            if(iflag[i][j][k+kk]==1){
              kappa_sum += kappa[i][j][k+kk];
              flag +=1.0;
            }
          }
          kappa[i][j][k] = kappa_sum / flag;
          iflagx[i][j][k] = 1;
        }
        sum  = 0;
        flag = 0.0;
        kappa_sum = 0.0; 
      } 
    }
    for_aijk(i,j,k){
      iflag[i][j][k] = iflagx[i][j][k];
    } 
  }
  for_aijk(i,j,k) {
    iflag[i][j][k]=stmp2[i][j][k];
    iflagx[i][j][k]=stmp3[i][j][k];
  }


  iflagx.exchange();
  iflag.exchange();
  
  return;
}
# endif

/* 1 to ej() 2D */

#if 1
void VOF::curv_HF_ext(){

  int sum=0;
  real flag=0.0;
  real kappa_sum=0.0;

  for_aijk(i,j,k) {
    stmp2[i][j][k]=iflag[i][j][k];
    stmp3[i][j][k]=iflagx[i][j][k];
    iflag[i][j][k] = iflagx[i][j][k];
  }
  iflag.exchange();

  if(time->current_step()==2){
    std::cout<<"iflag[10][3][16]= "<<iflag[10][3][16]<<"\n";
    std::cout<<"iflagx[10][3][16]= "<<iflagx[10][3][16]<<"\n";
  }

  for(int layer=1; layer<=3; layer++){
    for(int i=2; i<=ei()-1; i++){
      for(int j=1; j<=ej(); j++){
        for(int k=2; k<=ek()-1; k++){ 
          if(iflag[i][j][k]==0){
            for(int ii=-1; ii<=1; ii++){
              sum += iflag[i+ii][j][k];
            } 
            for(int kk=-1; kk<=1; kk++){
              sum += iflag[i][j][k+kk];
            }
#if 1
            if(time->current_step()==2 && i==10 && j==3 && k==17 && layer==1){
              std::cout<<"sum= "<<sum<<"\n";
              std::cout<<"iflag[9][3][17]= "<<iflag[9][3][17]<<"\n";
              std::cout<<"kappa[9][3][17]= "<<kappa[9][3][17]<<"\n";
              std::cout<<"iflag[10][3][17]= "<<iflag[10][3][17]<<"\n";
              std::cout<<"kappa[10][3][17]= "<<kappa[10][3][17]<<"\n";
              std::cout<<"iflag[11][3][17]= "<<iflag[11][3][17]<<"\n";
              std::cout<<"kappa[11][3][17]= "<<kappa[11][3][17]<<"\n";
              std::cout<<"iflag[10][3][16]= "<<iflag[10][3][16]<<"\n";
              std::cout<<"kappa[10][3][16]= "<<kappa[10][3][16]<<"\n";
              std::cout<<"iflag[10][3][18]= "<<iflag[10][3][18]<<"\n";
              std::cout<<"kappa[10][3][18]= "<<kappa[10][3][18]<<"\n";
            }
            if(time->current_step()==2 && i==10 && j==3 && k==16 && layer==1){
              std::cout<<"sum= "<<sum<<"\n";
              std::cout<<"iflag[9][3][16]= "<<iflag[9][3][16]<<"\n";
              std::cout<<"kappa[9][3][16]= "<<kappa[9][3][16]<<"\n";
              std::cout<<"iflag[10][3][16]= "<<iflag[10][3][16]<<"\n";
              std::cout<<"kappa[10][3][16]= "<<kappa[10][3][16]<<"\n";
              std::cout<<"iflag[11][3][16]= "<<iflag[11][3][16]<<"\n";
              std::cout<<"kappa[11][3][16]= "<<kappa[11][3][16]<<"\n";
              std::cout<<"iflag[10][3][15]= "<<iflag[10][3][15]<<"\n";
              std::cout<<"kappa[10][3][15]= "<<kappa[10][3][15]<<"\n";
              std::cout<<"iflag[10][3][17]= "<<iflag[10][3][17]<<"\n";
              std::cout<<"kappa[10][3][17]= "<<kappa[10][3][17]<<"\n";
            }

#endif

             
            if(sum>0){
              for(int ii=-1; ii<=1; ii++){
                if(iflag[i+ii][j][k]==1){
                  kappa_sum += kappa[i+ii][j][k];
                  flag += 1.0; 
                }
              } 
              for(int kk=-1; kk<=1; kk++){
                if(iflag[i][j][k+kk]==1){
                  kappa_sum += kappa[i][j][k+kk];
                  flag +=1.0;
                }
              }
#if 1
            if(time->current_step()==2 && i==10 && j==3 && k==16 && layer==1){
              std::cout<<"kappa_sum= "<<kappa_sum<<"flag= "<<flag<<"\n";
            }
#endif  
              kappa[i][j][k] = kappa_sum / flag;
              iflagx[i][j][k] = 1;
            }
            sum  = 0;
            flag = 0.0;
            kappa_sum = 0.0; 
          }
        }
      } 
    }
    for_aijk(i,j,k){
      iflag[i][j][k] = iflagx[i][j][k];
    }
    kappa.exchange();
    iflag.exchange(); 
  }
  for_aijk(i,j,k) {
    iflag[i][j][k]=stmp2[i][j][k];
    iflagx[i][j][k]=stmp3[i][j][k];
  }


  iflagx.exchange();
  iflag.exchange();
  
  return;
}
# endif



