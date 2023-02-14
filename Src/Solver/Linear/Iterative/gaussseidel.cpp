#include "gaussseidel.h"

//#define DEBUG

/***************************************************************************//**
*  \brief Implementation of the Gauss-Seidel solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
bool GaussSeidel :: solve(Matrix & A, Scalar & x, Scalar & b,
                     const MinIter & mini,
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
  boil::oout<<"GaussSeidel: iteration= 0  residual= "<<sqrt(res)<<" ratio= "
            <<sqrt(res/res0)<<" mini= "<<mini<<" maxi= "<<mi<<" res_tol= "
            <<res_tol<<"\n";

  /* should res be scaled with A and x? */
  if(sqrt(res) < res_tol) return true; // temporary meassure

  bool converged(false);
  int it;
  for(it=0; it<mi; it++) {
    
    /* perform one iteration step */
    for_vijk(x,i,j,k) {
      x[i][j][k] = A.ci[i][j][k]*(r[i][j][k]+A.c[i][j][k]*x[i][j][k]);
    }

    //std::cout<<"gaussseidel: "<<x[3][3][0+2]<<" "<<x[3][3][1+2]<<" "<<A.b[3][3][1+2]<<" "<<b[3][3][1+2]<<"\n";

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
