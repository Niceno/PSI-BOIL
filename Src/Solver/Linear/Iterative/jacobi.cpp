#include "jacobi.h"

//#define DEBUG

/***************************************************************************//**
*  \brief Implementation of the Jacobi solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
bool Jacobi :: solve(Matrix & A, Scalar & x, Scalar & b, const MinIter & mini,
                     const MaxIter & mi,
                     const char * name,
                     const ResRat & res_rat, const ResTol & res_tol,
                     const real scale,
                     const int stalecount,
                     const bool precform) {

  r = x.shape(); r=0.0;

  /*----------------------+
  |  compute r = b - A x  |
  +----------------------*/
  r = b - A * x;
  real res = r.dot(r); 
  real res0 = res;

#ifdef DEBUG
  OMS(------------);
  OPR(sqrt(res0));
  OPR(res_tol);
  OPR(res_rat);
#endif
  boil::oout<<"Jacobi: iteration= 0  residual= "<<sqrt(res)<<" ratio= "
            <<sqrt(res/res0)<<" mini= "<<mini<<" maxi= "<<mi<<" res_tol= "
            <<res_tol<<"\n";

  /* should res be scaled with A and x? */
  if(sqrt(res) < res_tol) return true; // temporary meassure

  bool converged(false);
  int it;
  for(it=0; it<mi; it++) {
    
    /* perform one iteration step */
    for_vijk(x,i,j,k) {
      // Lubomir: this must be Gauss-Seidel
      //x[i][j][k] = A.ci[i][j][k]*(r[i][j][k]+A.c[i][j][k]*x[i][j][k]);
      // Yohei
      x[i][j][k] = A.ci[i][j][k]*(
                       b[i][j][k]
                     + A.w[i][j][k] * x[i-1][j][k]
                     + A.e[i][j][k] * x[i+1][j][k]
                     + A.s[i][j][k] * x[i][j-1][k]
                     + A.n[i][j][k] * x[i][j+1][k]
                     + A.b[i][j][k] * x[i][j][k-1]
                     + A.t[i][j][k] * x[i][j][k+1]
                    );
    }
    //boil::oout<<"jacobi: "<<x[3][3][0+2]<<" "<<x[3][3][1+2]<<" "<<A.c[3][3][1+2]<<" "<<A.ci[3][3][3]<<" "<<A.t[3][3][1+2]<<" "<<b[3][3][1+2]<<"\n";

    /* exchange and compute residual */
    r = b - A * x;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    res = r.dot(r);

#ifdef DEBUG
    OPR( sqrt(res) );
#endif
    if(it >= 30 && it%10 ==0)
      boil::oout<<"iteration= "<<it<<" residual= "<<sqrt(res)
                <<" ratio= "<<sqrt(res/res0)<<"\n";

    /* should res be scaled with A and x? */
    if( sqrt(res) < res_tol && it >= mini-1 ) { converged = true; break; }

    if( sqrt(res) < sqrt(res0) * res_rat && it >= mini-1 ) { converged = true; break; } 
  }

  /* for normalisation */
  r = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << sqrt(res/r.dot(r)) 
                            << ", ratio = " << sqrt ( res/res0 )
                            << ", iterations = " << it+1 
                            << boil::endl;

  return converged;

}
