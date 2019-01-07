#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_cc(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Centered Columns method in E.Aulisa,JCP,225(2007),2301-2319
*         Resluts: nx, ny, nz
*******************************************************************************/
#if 0
  for_ijk(i,j,k) {
    real mxX, myX, mzX;
    mxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
    myX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
                - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k])); 
    mzX = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
                - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 
    normalize(mxX,myX,mzX);

#if 1
    if(i==13 && j==21 && k==25){
      std::cout<<"mxX="<<mxX<<" myX="<<myX<<" mzX="<<mzX<<"\n";
    }
#endif

    real mxY, myY, mzY;
    mxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
                - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
    myY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
    mzY = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));
    normalize(mxY,myY,mzY);

#if 1
    if(i==13 && j==21 && k==25){
      std::cout<<"mxY="<<mxY<<" myY="<<myY<<" mzY="<<mzY<<"\n";
    }
#endif

    real mxZ, myZ, mzZ;
    mxZ = 0.5 * ( (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
                - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
    myZ = 0.5 * ( (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
    mzZ = copysign(1.0,+(sca[i][j][k+1]-sca[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

#if 1
    if(i==13 && j==21 && k==25){
      std::cout<<"mxZ="<<mxZ<<" myZ="<<myZ<<" mzZ="<<mzZ<<"\n";
    }
#endif

    if (fabs(mxX)<fabs(myY)) {
      if (fabs(myY)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxY;
        ny[i][j][k]=myY;
        nz[i][j][k]=mzY;
      }
    } else {
      if (fabs(mxX)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxX;
        ny[i][j][k]=myX;
        nz[i][j][k]=mzX;
      }
    }
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  //insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
  insert_bc_norm_cc(sca);

  /* normalize */
  //for_avijk(sca,i,j,k) {
  //  normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  //}

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

#if 1
  if(time->current_step() % 858 == 0) {
    boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
    //exit(0);
  }
#endif
  return;
#endif

for_ijk(i,j,k) {
    real mxX, myX, mzX, myXb, myXc, myXf, mzXb, mzXc, mzXf;
    mxX = copysign(1.0,+(sca[i+1][j-1][k-1]+sca[i+1][j-1][k]+sca[i+1][j-1][k+1]+
                         sca[i+1][j][k-1]  +sca[i+1][j][k]  +sca[i+1][j][k+1]+
                         sca[i+1][j+1][k-1]+sca[i+1][j+1][k]+sca[i+1][j+1][k+1]- 
                         sca[i-1][j-1][k-1]-sca[i-1][j-1][k]-sca[i-1][j-1][k+1]-
                         sca[i-1][j][k-1]  -sca[i-1][j][k]  -sca[i-1][j][k+1]-  
                         sca[i-1][j+1][k-1]-sca[i-1][j+1][k]-sca[i-1][j+1][k+1]));
    myXb = (sca[i+1][j][k]   +sca[i][j][k]  +sca[i-1][j][k])
          -(sca[i+1][j-1][k] +sca[i][j-1][k]+sca[i-1][j-1][k]); 
    myXc = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
                 - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k]));
    myXf = (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k]) 
          -(sca[i+1][j][k]  +sca[i][j][k]  +sca[i-1][j][k]);
    mzXb = (sca[i+1][j][k]  +sca[i][j][k]  +sca[i-1][j][k])
          -(sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1]);
    mzXc = 0.5 * ( (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
                 - (sca[i+1][j][k-1]+sca[i][j][k-1]+sca[i-1][j][k-1])); 
    mzXf = (sca[i+1][j][k+1]+sca[i][j][k+1]+sca[i-1][j][k+1])
          -(sca[i+1][j][k]  +sca[i][j][k]  +sca[i-1][j][k]);

    if (fabs(myXb)<fabs(myXc)) {
      if (fabs(myXc)<fabs(myXf)){
        myX = myXf;
      } else{
        myX = myXc;
      }
    } else{
      if(fabs(myXb)<fabs(myXf)) {
        myX = myXf;
      } else{
        myX = myXb;
      }
    }

    if (fabs(mzXb)<fabs(mzXc)) {
      if (fabs(mzXc)<fabs(mzXf)){
        mzX = mzXf;
      } else{
        mzX = mzXc;
      }
    } else{
      if(fabs(mzXb)<fabs(mzXf)) {
        mzX = mzXf;
      } else{
        mzX = mzXb;
      }
    }

    normalize(mxX,myX,mzX);

#if 0
    if(i==13 && j==21 && k==25){
      std::cout<<"mxX="<<mxX<<" myX="<<myX<<" mzX="<<mzX<<"\n";
    }
#endif

    real mxY, myY, mzY, mxYb, mxYc, mxYf, mzYb, mzYc, mzYf;
    mxYb = (sca[i][j-1][k]  +sca[i][j][k]  +sca[i][j+1][k])
          -(sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k]); 
    mxYc = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
                 - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k])); 
    mxYf = (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
          -(sca[i][j-1][k]  +sca[i][j][k]  +sca[i][j+1][k]);
    myY = copysign(1.0,+(sca[i-1][j+1][k-1]+sca[i-1][j+1][k]+sca[i-1][j+1][k+1]+
                         sca[i][j+1][k-1]  +sca[i][j+1][k]  +sca[i][j+1][k+1]+
                         sca[i+1][j+1][k-1]+sca[i+1][j+1][k]+sca[i+1][j+1][k+1]- 
                         sca[i-1][j-1][k-1]-sca[i-1][j-1][k]-sca[i-1][j-1][k+1]-
                         sca[i][j-1][k-1]  -sca[i][j-1][k]  -sca[i][j-1][k+1]-
                         sca[i+1][j-1][k-1]-sca[i+1][j-1][k]-sca[i+1][j-1][k+1]));
    mzYb = (sca[i][j-1][k]  +sca[i][j][k]  +sca[i][j+1][k])
          -(sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]);
    mzYc = 0.5 * ( (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j][k-1]+sca[i][j+1][k-1]));
    mzYf = (sca[i][j-1][k+1]+sca[i][j][k+1]+sca[i][j+1][k+1])
          -(sca[i][j-1][k]  +sca[i][j][k]  +sca[i][j+1][k]);

    if (fabs(mxYb)<fabs(mxYc)) {
      if (fabs(mxYc)<fabs(mxYf)){
        mxY = mxYf;
      } else{
        mxY = mxYc;
      }
    } else{
      if(fabs(mxYb)<fabs(mxYf)) {
        mxY = mxYf;
      } else{
        mxY = mxYb;
      }
    }

    if (fabs(mzYb)<fabs(mzYc)) {
      if (fabs(mzYc)<fabs(mzYf)){
        mzY = mzYf;
      } else{
        mzY = mzYc;
      }
    } else{
      if(fabs(mzYb)<fabs(mzYf)) {
        mzY = mzYf;
      } else{
        mzY = mzYb;
      }
    }

    normalize(mxY,myY,mzY);

#if 0
    if(i==13 && j==21 && k==25){
      std::cout<<"mxY="<<mxY<<" myY="<<myY<<" mzY="<<mzY<<"\n";
    }
#endif

    real mxZ, myZ, mzZ, mxZb, mxZc, mxZf, myZb, myZc, myZf;
    mxZb = (sca[i][j][k-1]  +sca[i][j][k]  +sca[i][j][k+1])
          -(sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1]);
    mxZc = 0.5 * ((sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
                - (sca[i-1][j][k-1]+sca[i-1][j][k]+sca[i-1][j][k+1])); 
    mxZf = (sca[i+1][j][k-1]+sca[i+1][j][k]+sca[i+1][j][k+1])
          -(sca[i][j][k-1]  +sca[i][j][k]  +sca[i][j][k+1]);
    myZb = (sca[i][j][k-1]  +sca[i][j][k]  +sca[i][j][k+1])
          -(sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1]);
    myZc = 0.5 * ((sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
                - (sca[i][j-1][k-1]+sca[i][j-1][k]+sca[i][j-1][k+1])); 
    myZf = (sca[i][j+1][k-1]+sca[i][j+1][k]+sca[i][j+1][k+1])
          -(sca[i][j][k-1]  +sca[i][j][k]  +sca[i][j][k+1]);
    mzZ = copysign(1.0,+(sca[i-1][j-1][k+1]+sca[i-1][j][k+1]+sca[i-1][j+1][k+1]+
                         sca[i][j-1][k+1]  +sca[i][j][k+1]  +sca[i][j+1][k+1]+
                         sca[i+1][j-1][k+1]+sca[i+1][j][k+1]+sca[i+1][j+1][k+1]-
                         sca[i-1][j-1][k-1]-sca[i-1][j][k-1]-sca[i-1][j+1][k-1]-
                         sca[i][j-1][k-1]  -sca[i][j][k-1]  -sca[i][j+1][k-1]-
                         sca[i+1][j-1][k-1]-sca[i+1][j][k-1]-sca[i+1][j+1][k-1]));
  
    if (fabs(mxZb)<fabs(mxZc)) {
      if (fabs(mxZc)<fabs(mxZf)){
        mxZ = mxZf;
      } else{
        mxZ = mxZc;
      }
    } else{
      if(fabs(mxZb)<fabs(mxZf)) {
        mxZ = mxZf;
      } else{
        mxZ = mxZb;
      }
    }

    if (fabs(myZb)<fabs(myZc)) {
      if (fabs(myZc)<fabs(myZf)){
        myZ = myZf;
      } else{
        myZ = myZc;
      }
    } else{
      if(fabs(myZb)<fabs(myZf)) {
        myZ = myZf;
      } else{
        myZ = myZb;
      }
    }


    normalize(mxZ,myZ,mzZ);

#if 0
    if(i==13 && j==21 && k==25){
      std::cout<<"mxZ="<<mxZ<<" myZ="<<myZ<<" mzZ="<<mzZ<<"\n";
    }
#endif

    if (fabs(mxX)<fabs(myY)) {
      if (fabs(myY)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxY;
        ny[i][j][k]=myY;
        nz[i][j][k]=mzY;
      }
    } else {
      if (fabs(mxX)<fabs(mzZ)) {
        nx[i][j][k]=mxZ;
        ny[i][j][k]=myZ;
        nz[i][j][k]=mzZ;
      } else {
        nx[i][j][k]=mxX;
        ny[i][j][k]=myX;
        nz[i][j][k]=mzX;
      }
    }
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  //insert_bc_gradphic(sca); 
  
  /* normal vector on boundary plane */
    insert_bc_norm_cc(sca);
  
  //      /* normalize */
  // for_avijk(sca,i,j,k) {
  //  normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
   //}
  
    nx.exchange_all();
    ny.exchange_all();
    nz.exchange_all();
  
  #if 0
    if(time->current_step() % 858 == 0) {
      boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
      //exit(0);
    }
  #endif
    return;
}
