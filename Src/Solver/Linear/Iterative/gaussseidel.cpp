#include "gaussseidel.h"

//#define DEBUG

/***************************************************************************//**
*  \brief Implementation of the Gauss-Seidel solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
void GaussSeidel :: solve(Matrix & A, Scalar & x, Scalar & b,
                          const MinIter & minit, const MaxIter & maxit,
                          const char * name,
                          const ResRat & res_rat, const ResTol & res_tol) {
                          //const real scale,
                          //const int stalecount,
                          //const bool precform) {

  r = x.shape(); r=0.0;

  /*-----------------------------------+
  |  exchange and compute r = b - A x  |
  +-----------------------------------*/
  r = b - A * x;
  //real res = sqrt(r.dot_voldiv_avg(r))/scale;
  real res = r.dot(r);
  real res0 = res;

  /* staleness vector */
  //std::vector<real> resvect;
  //if(stalecount>0) {
  //  resvect.resize(stalecount);
  //  for(auto & r : resvect)
  //    r = boil::unreal;
  //}

#ifdef DEBUG
  OMS(------------);
  OPR(res0);
  OPR(res_tol);
  OPR(res_rat);
#endif

  /* should res be scaled with A and x? */
  if(res < res_tol) return; // temporary meassure

  //bool converged(false);
  int it;
  for(it=0; it<maxit; it++) {
    
    /* perform one iteration step */
    for_vijk(x,i,j,k) {
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

    /* exchange and compute residual */
    r = b - A * x;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    //res = sqrt(r.dot_voldiv_avg(r))/scale;
    res = r.dot(r);

#ifdef DEBUG
    OPR( res );
    if(stalecount>0) {
      boil::oout<<i<<" ";
      for(auto & r : resvect)
        boil::oout<<" "<<r;
      boil::oout<<" "<<res<<boil::endl;
    }
#endif

    /* should res be scaled with A and x? */
    if( res < res_tol && it >= minit-1 )  break;

    if( res < res0 * res_rat && it >= minit-1 ) break;

#if 0
    if(stalecount>0) {
      bool staleflag(true);
      for(auto & r : resvect) {
        if(res<r) {
          staleflag = false;
          break;
        }
      }
      if(staleflag) {
        if(name!=NULL) {
          boil::oout << name  << " staled!";
          for(auto & r : resvect)
            boil::oout<<" "<<r;
          boil::oout<<" "<<res<<boil::endl;
        }
        break;
      } else {
        std::rotate(resvect.begin(),resvect.begin()+1,resvect.end());
        resvect.back() = res;
      }
    }
#endif
  }

  x.exchange();

  /* for normalisation */
  r = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << sqrt(res/r.dot(r))
                            << ", ratio = " << sqrt ( res/res0 )
                            << ", iterations = " << it+1 
                            << boil::endl;

  return;

}
