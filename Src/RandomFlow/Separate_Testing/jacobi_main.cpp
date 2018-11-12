/* compile with: g++ jacobi_main.cpp ../Global/global_malloc.cpp */
#include <cmath>
#include <iostream>

#include "../Global/global_malloc.h"

namespace {
  inline void rot(real ** a, const real s, const real tau, const int i,
    const int j, const int k, const int l)
  {
    real g,h;

    g=a[i][j];
    h=a[k][l];
    a[i][j]=g-s*(h+g*tau);
    a[k][l]=h+s*(g-h*tau);
  }
}

/******************************************************************************/
void jacobi(real ** a, real * d, real ** v, const int n, int * nrot)
{
  int i,j,ip,iq;
  real tresh,theta,tau,t,sm,s,h,g,c;

  real * b = new real[n];
  real * z = new real[n];
  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0)
      return;
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
          && (fabs(d[iq])+g) == fabs(d[iq]))
            a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((fabs(h)+g) == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0;j<ip;j++)
            rot(a,s,tau,j,ip,j,iq);
          for (j=ip+1;j<iq;j++)
            rot(a,s,tau,ip,j,j,iq);
          for (j=iq+1;j<n;j++)
            rot(a,s,tau,ip,j,iq,j);
          for (j=0;j<n;j++)
            rot(v,s,tau,j,ip,j,iq);
          ++(*nrot);
        }
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
 
  delete [] b;
  delete [] z;

  std::cout << "Too many iterations in routine jacobi" << std::endl;
}

/******************************************************************************/
main() {

  const int N=3;
  int nrot;

  real ** A;
 
  real * B, ** V;

  alloc2d(&A,N,N);
  alloc1d(&B,N);
  alloc2d(&V,N,N);

/*
  A[0][0] = 1.0;   A[0][1] = A[1][0] = -0.5;   A[0][2] = A[2][0] = 0.0;
  A[1][1] = 1.0;   A[1][2] = A[2][1] = -0.5;
  A[2][2] = 1.0;
*/

  A[0][0] = 1.0;   A[0][1] = A[1][0] =  0.0;   A[0][2] = A[2][0] = 0.0;
  A[1][1] = 0.2;   A[1][2] = A[2][1] =  0.0;
  A[2][2] = 0.3;

  for(int i=0; i<N; i++) 
    std::cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << std::endl;

  for(int i=0; i<N; i++) 
    std::cout << V[i][0] << " " << V[i][1] << " " << V[i][2] << std::endl;

  jacobi(A, B, V, N, & nrot);

  for(int i=0; i<N; i++) 
    std::cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << std::endl;

  for(int i=0; i<N; i++) 
    std::cout << V[i][0] << " " << V[i][1] << " " << V[i][2] << std::endl;

  for(int i=0; i<N; i++) 
    std::cout << B[i] << std::endl;

  std::cout << "nrot = " << nrot << std::endl;

}
