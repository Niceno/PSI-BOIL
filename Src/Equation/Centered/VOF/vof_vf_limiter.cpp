#include "vof.h"
/******************************************************************************/
/*  If volume fraction is larger than one or smaller than zero, spread the    */
/*  unwanted volume fraction to the neighbour cells with volume fraction      */
/*  between  zero and one.                                                    */
/******************************************************************************/
#if 0
void VOF::vf_limiter(){

  int count=0;
  for_ijk(i,j,k){
    if(phi[i][j][k]>1.0){
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              count++;
            }
          }
        }
      }
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              phi[i+ii][j+jj][k+kk] += (phi[i][j][k]-1.0)/count;
            }
          }
        }
      }
      std::cout.setf(std::ios_base::scientific);
      std::cout<<"(phi[i][j][k]-1.0)/count="<<(phi[i][j][k]-1.0)/count<<"\n";
      std::cout.unsetf(std::ios_base::floatfield);
      if(count>0){
        phi[i][j][k] = 1.0;
      }
      count = 0;
    }
    if(phi[i][j][k]<0.0){
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              count++;
            } 
          }
        }
      }
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              phi[i+ii][j+jj][k+kk] -= (0.0-phi[i][j][k])/count;
            }
          }
        }
      }
      if(count>0){
        phi[i][j][k] = 0.0;
      }
      count = 0;
    }
  }
  return;
}
#endif

#if 0
void VOF::vf_limiter(){
  int a,b,c;
  int count=0;
  for_ijk(i,j,k){
    if(phi[i][j][k]>1.0){
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              a = i+ii;
              b = j+jj;
              c = k+kk;
              count = 1;
            }
          }
        }
      }
      if(count=1){
        phi[a][b][c] += phi[i][j][k]-1.0;
        phi[i][j][k]  = 1.0;
        a = 0;
        b = 0;
        c = 0;
        count = 0;
      }
    }

    if(phi[i][j][k]<0.0){
      for(int ii=-1; ii<=1; ii++){
        for(int jj=-1; jj<=1; jj++){
          for(int kk=-1; kk<=1; kk++){
            if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
              a = i+ii;
              b = j+jj;
              c = k+kk;
              count = 1;
            } 
          }
        }
      }
      if(count=1){
        phi[a][b][c] -= 0.0-phi[i][j][k];
        phi[i][j][k]  = 0.0;
        a = 0;
        b = 0;
        c = 0;
        count = 0;
      }
    }
  }
  return;
}
#endif

#if 0
void VOF::vf_limiter(){

  int count=0;
  for(int i=2; i<=ei()-1; i++){  
    for(int j=2; j<=ej()-1; j++){
      for(int k=2; k<=ek()-1; k++){
        if(phi[i][j][k]>1.0){
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  count++;
                }
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  phi[i+ii][j+jj][k+kk] += (phi[i][j][k]-1.0)/count;
                }
              }
            }
          }
          std::cout.setf(std::ios_base::scientific);
          std::cout<<"(phi[i][j][k]-1.0)/count="<<(phi[i][j][k]-1.0)/count<<"\n";
          std::cout.unsetf(std::ios_base::floatfield);
          if(count>0){
            phi[i][j][k] = 1.0;
          }
          count = 0;
        }
        if(phi[i][j][k]<0.0){
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  count++;
                }
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  phi[i+ii][j+jj][k+kk] -= (0.0-phi[i][j][k])/count;
                }
              }
            }
          }
          if(count>0){
            phi[i][j][k] = 0.0;
          }
          count = 0;
        }
      }
    }
  }
  return;
}
#endif

//1 to ej()


#if 0
void VOF::vf_limiter(){

  int count=0;
  for(int i=2; i<=ei()-1; i++){  
    for(int j=1; j<=ej(); j++){
      for(int k=2; k<=ek()-1; k++){
        if(phi[i][j][k]>1.0){
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  count++;
                }
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  phi[i+ii][j+jj][k+kk] += (phi[i][j][k]-1.0)/count;
                }
              }
            }
          }
          std::cout.setf(std::ios_base::scientific);
          std::cout<<"(phi[i][j][k]-1.0)/count="<<(phi[i][j][k]-1.0)/count<<"\n";
          std::cout.unsetf(std::ios_base::floatfield);
          if(count>0){
            phi[i][j][k] = 1.0;
          }
          count = 0;
        }
        if(phi[i][j][k]<0.0){
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  count++;
                }
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int jj=-1; jj<=1; jj++){
              for(int kk=-1; kk<=1; kk++){
                if(boil::micro<phi[i+ii][j+jj][k+kk] && phi[i+ii][j+jj][k+kk]<(1.0-boil::micro)){
                  phi[i+ii][j+jj][k+kk] -= (0.0-phi[i][j][k])/count;
                }
              }
            }
          }
          if(count>0){
            phi[i][j][k] = 0.0;
          }
          count = 0;
        }
      }
    }
  }
  return;
}
#endif

/* 1 to ej() 2D */

//The start and end of the loop need to be reconsidered if taking into account boundary condition.
//This only works for 2D. 


#if 1
void VOF::vf_limiter(){

  int count=0;
  for(int i=2; i<=ei()-1; i++){  
    for(int j=1; j<=ej(); j++){
      for(int k=2; k<=ek()-1; k++){
        if(phi[i][j][k]>1.0){
          for(int ii=-1; ii<=1; ii++){
            for(int kk=-1; kk<=1; kk++){
              if(boil::micro<phi[i+ii][j][k+kk] && phi[i+ii][j][k+kk]<(1.0-boil::micro)){
                count++;
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int kk=-1; kk<=1; kk++){
              if(boil::micro<phi[i+ii][j][k+kk] && phi[i+ii][j][k+kk]<(1.0-boil::micro)){
                phi[i+ii][j][k+kk] += (phi[i][j][k]-1.0)/count;
              }
            }
          }
#if 0
          std::cout.setf(std::ios_base::scientific);
          std::cout<<"(phi[i][j][k]-1.0)/count="<<(phi[i][j][k]-1.0)/count<<"\n";
          std::cout.unsetf(std::ios_base::floatfield);
#endif
          if(count>0){
            phi[i][j][k] = 1.0;
          }
          count = 0;
        }
        if(phi[i][j][k]<0.0){
          for(int ii=-1; ii<=1; ii++){
            for(int kk=-1; kk<=1; kk++){
              if(boil::micro<phi[i+ii][j][k+kk] && phi[i+ii][j][k+kk]<(1.0-boil::micro)){
                count++;
              }
            }
          }
          for(int ii=-1; ii<=1; ii++){
            for(int kk=-1; kk<=1; kk++){
              if(boil::micro<phi[i+ii][j][k+kk] && phi[i+ii][j][k+kk]<(1.0-boil::micro)){
                phi[i+ii][j][k+kk] -= (0.0-phi[i][j][k])/count;
              }
            }
          }
          if(count>0){
            phi[i][j][k] = 0.0;
          }
          count = 0;
        }
      }
    }
  }
  return;
}
#endif

