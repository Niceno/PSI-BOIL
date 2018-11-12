#include "randomflow.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>

/******************************************************************************/
real RandomFlow::sclp(const real * A, const real * B) {
  return (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]);
}

/******************************************************************************/
void RandomFlow::vecp(real * A, const real * B, const real * C) {
  A[0]=B[1]*C[2]-B[2]*C[1];
  A[1]=B[2]*C[0]-B[0]*C[2];
  A[2]=B[0]*C[1]-B[1]*C[0];
}

/******************************************************************************/
real RandomFlow::gauss()
/*--------------------------------------+
|  From Numerical Recepes in C, Ch.7.2  |
|  Zero mean and unit variance          |
+--------------------------------------*/
{
  static int iset = 0;
  static real gset;
  real fac,rsq,v1,v2;
  if (iset == 0)
  {  
    do
    {
      v1=2.*boil::random_number()-1.;
      v2=2.*boil::random_number()-1.;
      rsq=v1*v1+v2*v2;
    }  while (rsq >= 1.0 || rsq == 0);
    fac=(real)sqrt(-2.*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }
  else
  {
    iset=0;
    return gset;
  }
}

/******************************************************************************/
void RandomFlow::gaussn(real **Y, real d, int n, int m) {
/*-------------------------------------------------+
|  Generates an array of n random numbers with a   |
|  Normal (Gaussian) distribution with a standard  |
|  deviation of d, and stores them in array Y.     |
+-------------------------------------------------*/
  for(int i=0; i<n; i++) 
    for(int j=0; j<m; j++) 
      Y[i][j] = d * gauss();
}

/******************************************************************************/
void RandomFlow::gaussn(real *Y, real d, int n) {
/*-------------------------------------------------+
|  Generates an array of n random numbers with a   |
|  Normal (Gaussian) distribution with a standard  |
|  deviation of d, and stores them in array Y.     |
+-------------------------------------------------*/
  real *x;
  for (x = Y; x - Y < n; x++)
    *x = d*gauss();
}
