#include <cmath>
#include <cstdio>
#include <vector>

#include "definitions.h"
#include "variables.h"

/******************************************************************************/
int solve_cg(const matrix            & A, 
                   std::vector<real> & x, 
             const std::vector<real> & b, 
             const real                tol, 
             const int                 prec, 
             const int                 niter) {
/*---------------------------------------------------------------------+
|   this routine solves the linear systems of equations [A]{x}={b}     |
|   by the conjugate gradient method.                                  |
|                                                                      |
|   this routine allows to precondition the system of equations        |
|   with:                                                              |
|     1. diagonal preconditioning                                      |
|     2. incomplete cholesky preconditioning                           |
|                                                                      |
|   the type of precondtioning is chosen by setting the variable prec  |
|   to 0 (no preconditioning), 1 (diagonal preconditioning) or 2       |
|   (incomplete cholesky preconditioning)                              |
+---------------------------------------------------------------------*/
  real              rho, rho_old;
  std::vector<real> Ax, D, p, q, r, z;
 
  const int N = A.size();

  Ax.resize(N);
  D .resize(N);
  p .resize(N);
  q .resize(N);
  r .resize(N);
  z .resize(N);

  /*------------------+
  |  preconditioning  |
  +------------------*/

  /* 1) diagonal preconditioning */
  if(prec==1)         
    for(int i=0;i<N;i++)
      D[i]=A.val[i][0];

  /* 2) incomplete cholesky preconditioning */
  else if(prec==2)
    for(int i=0;i<N;i++) {
      real sum = A.val[i][0];
      for(int j=1; j<=A.con[i][0]; j++) {
        const int k=A.con[i][j];
        if(k<i) {    /* lower triangular matrix */
          sum = sum - D[k] * A.val[i][j]*A.val[i][j];               
        }
        D[i] = 1.0 / sum;
      }
    }

  /* .) no preconditioning */
  else                          
    for(int i=0;i<N;i++) {
      D[i]=1.0;
    }

  real bnrm2 = 0.0;
  for(int i=0;i<N;i++) {
    bnrm2 = bnrm2+b[i]*b[i];
  }  
  bnrm2 = sqrt(bnrm2);
  if(bnrm2 == 0.0) 
    bnrm2 = 1.0;  

  /*-------------+
  |  r = b - Ax  |
  +-------------*/
  for(int i=0;i<N;i++) {
    Ax[i]=A.val[i][0]*x[i];
    for(int j=1; j<=A.con[i][0]; j++) { 
      Ax[i] = Ax[i] + A.val[i][j] * x[A.con[i][j]];
    }
    r[i]=b[i]-Ax[i]; 
  }

  /*-----------------------------+
  |  calculate initial residual  |
  +-----------------------------*/
  real error=0.0;
  for(int i=0;i<N;i++) {
    error=error + r[i]*r[i];
  }
  error=sqrt(error)/bnrm2;
  if(error<tol) 
    return 0;

  /*----->------+
  |             |
  ^  main loop  v
  |             |
  +------<-----*/
  for(int iter=0; iter<niter; iter++) {
           
    /*----------------+  
    |  solve Mz  = r  |
    +----------------*/

    /* 1) diagonal preconditioning */   
    if(prec==1)
      for(int i=0;i<N;i++) {
        z[i]  = r[i]  / D[i];
      }

    /* 2) incomplete cholesky preconditioning */   
    else if(prec==2) {

      /* forward substitution */
      for(int i=0;i<N;i++) {
        real sum  = r[i];
        for(int j=1; j<=A.con[i][0]; j++) {
          const int k=A.con[i][j];
          if(k<i) { /* lower triangular matrix */
            sum  = sum  - A.val[i][j]*z[k];              
          }
        }
       z[i]  = sum  * D[i]; 
      }

      for(int i=0;i<N;i++) {
        z[i]  = z[i]  / ( D[i] + boil::atto );
      }

      /* backward substitution */
      for(int i=N-1;i>=0;i--) {
        real sum  = z[i];
        for(int j=1; j<=A.con[i][0]; j++) {
          const int k=A.con[i][j];
          if(k>i) { /* upper triangular matrix */ 
            sum  = sum  - A.val[i][j] * z[k];               
          }
        }
       z[i]  = sum  * D[i];
      }
    }

    /* .) no preconditioning */
    else
      for(int i=0; i<N; i++) 
	z[i]=r[i];

    /*--------------+
    |  rho = (z,r)  |
    +--------------*/
    rho = 0.0;
    for(int i=0;i<N;i++)
      rho = rho+z[i]*r[i];

    if(iter == 0) 
      for(int i=0;i<N;i++)
        p[i] = z[i]; 

    else {
      real beta = rho / rho_old;
      for(int i=0; i<N; i++) 
        p[i] = z[i] + beta * p[i];  
    }

    /*---------------------+
    |  q = A p             |
    |  alfa = (z,r)/(p,q)  |
    +---------------------*/
    real alfa = 0.0;
    for(int i=0;i<N;i++) {
      q[i]  = A.val[i][0] * p[i];
      for(int j=1; j<=A.con[i][0]; j++) { 
        q[i]  = q[i]  + A.val[i][j] * p[A.con[i][j]];
      }
      alfa = alfa + p[i] * q[i];
    }
    alfa = rho / alfa;

    /*-----------------+
    |  x = x + alfa p  |
    |  r = r - alfa q  |
    +-----------------*/
    for(int i=0;i<N;i++) {
      x[i]  = x[i]  + alfa*p[i];
      r[i]  = r[i]  - alfa*q[i];
    }
 
    /*--------------------+
    |  check convergence  |
    +--------------------*/
    error = 0.0;
    for(int i=0;i<N;i++)
      error=error+r[i]*r[i];

    error=sqrt(error) / bnrm2;
    if(error < tol) 
      return iter; 

    rho_old = rho;

  } /* iter */

 return niter;
}
