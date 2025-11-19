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

  /* should res be scaled with A and x? */
  if(sqrt(res) < res_tol) return true; // temporary meassure

  bool converged(false);
  int i;
  for(i=0; i<mi; i++) {
    
    /* perform one iteration step */
    for_vijk(x,i,j,k) {
      x[i][j][k] = A.ci[i][j][k]*(r[i][j][k]+A.c[i][j][k]*x[i][j][k]);
    }

    /* exchange and compute residual */
    r = b - A * x;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    res = r.dot(r);

#ifdef DEBUG
    OPR( sqrt(res) );
#endif

    /* should res be scaled with A and x? */
    if( sqrt(res) < res_tol && i >= mini-1 ) { converged = true; break; }

    if( sqrt(res) < sqrt(res0) * res_rat && i >= mini-1 ) { converged = true; break; } 
  }

  /* for normalisation */
  r = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << sqrt(res/r.dot(r)) 
                            << ", ratio = " << sqrt ( res/res0 )
                            << ", iterations = " << i+1 
                            << boil::endl;

  return converged;

}
