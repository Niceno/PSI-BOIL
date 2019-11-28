#include "vof.h"

//#define DEBUG

/******************************************************************************/
bool VOF::ev_solve(const ScalarInt & pflag, const Matrix & A,
                   const Scalar & b, Scalar & x, Scalar & xold, 
                   const bool init_guess, const int niter,
                   const ResRat & resrat) { 
/***************************************************************************//**
*  \brief Solve the problem A*x = b using Jacobi iterations.
*  pflag = -2,-1,1,-2 : solution domain
*******************************************************************************/

  /* initialize x */
  x = 0.;
  if(!init_guess) {
    xold = x;
  }

  /* Jacobi iterations */
  bool converged(false);
  real l2err_first, linferr_first;
  for(int n(0); n<niter; ++n) {
    for_ijk(i,j,k) {
      if(abs(pflag[i][j][k])<3) {
        x[i][j][k] = A.ci[i][j][k] 
                   * (b[i][j][k] - A.w[i][j][k]*xold[i-1][j][k]
                                 - A.e[i][j][k]*xold[i+1][j][k]
                                 - A.s[i][j][k]*xold[i][j-1][k]
                                 - A.n[i][j][k]*xold[i][j+1][k]
                                 - A.b[i][j][k]*xold[i][j][k-1]
                                 - A.t[i][j][k]*xold[i][j][k+1]);
      }
    }
    x.bnd_update();
    x.exchange();

    /* measure error */
    int nele(0);
    real l2err(0.0), linferr(0.0);
    for_ijk(i,j,k) {
      if(abs(pflag[i][j][k])<3) {
        nele++;
        real diff = fabs(x[i][j][k] - xold[i][j][k]);
        l2err += diff*diff;
        if(diff>linferr)
          linferr = diff;
      }
    }
    boil::cart.sum_int(&nele);
    boil::cart.sum_real(&l2err);
    boil::cart.max_real(&linferr);
    if(nele>0)
      l2err /= real(nele);
    l2err = sqrt(l2err);
    if(n==0) {
      l2err_first = l2err;
      linferr_first = linferr;
    }
#ifdef DEBUG
    boil::oout<<"VOF::ev_solve: "
              <<" "<<n<<" "<<linferr<<" "<<linferr/linferr_first<<" ( "<<nele<<" els)"
              <<boil::endl;
#endif
    xold = x;

    /* converged? */
    if(linferr/linferr_first<resrat) {
      converged = true;
#ifdef DEBUG
      boil::oout<<"VOF::ev_solve converged after "<<n<<" steps, final rel. error: "
                <<linferr/linferr_first<<"\n";
#endif
      break;
    }
  }

  if(!converged)
    boil::oout<<"VOF::ev_solve did not converge after "<<niter<<" steps! \n";

  return converged;
}
