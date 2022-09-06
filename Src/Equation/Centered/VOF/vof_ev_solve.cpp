#include "vof.h"

//#define DEBUG

/******************************************************************************/
bool VOF::ev_solve(const ScalarInt & pflag, const Matrix & A,
                   const Scalar & b, Scalar & x, Scalar & xold, 
                   const bool init_guess, const int niter,
                   const ResTol & restol) { 
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
  //real l2err_first, l2err(0.0);
  real linferr(0.0), linferr_first;
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
    //l2err=0.;
    linferr=0.;
    for_ijk(i,j,k) {
#if 0
      if(abs(pflag[i][j][k])<3) {
        nele++;
        real diff = fabs(x[i][j][k] - xold[i][j][k]);
        l2err += diff*diff;
        if(diff>linferr)
          linferr = diff;
      }
#elif 0
      if(abs(pflag[i][j][k])==1&&fabs(b[i][j][k])>boil::atto) {
        nele++;
        real diff = b[i][j][k] - A.w[i][j][k]*x[i-1][j][k]
                               - A.e[i][j][k]*x[i+1][j][k]
                               - A.s[i][j][k]*x[i][j-1][k]
                               - A.n[i][j][k]*x[i][j+1][k]
                               - A.b[i][j][k]*x[i][j][k-1]
                               - A.t[i][j][k]*x[i][j][k+1]
                               - A.c[i][j][k]*x[i][j][k];
        diff /= b[i][j][k];
        if(fabs(diff)>linferr)
          linferr = fabs(diff);
      }
#else
      if(abs(pflag[i][j][k])==1) {
        nele++;
        real urat = b[i][j][k] - A.w[i][j][k]*x[i-1][j][k]
                               - A.e[i][j][k]*x[i+1][j][k]
                               - A.s[i][j][k]*x[i][j-1][k]
                               - A.n[i][j][k]*x[i][j+1][k]
                               - A.b[i][j][k]*x[i][j][k-1]
                               - A.t[i][j][k]*x[i][j][k+1]
                               - A.c[i][j][k]*x[i][j][k];

        urat = fabs(urat)/dV(i,j,k)*time->dt();
        if(urat>linferr)
          linferr = urat;
      }
#endif
    }
    boil::cart.sum_int(&nele);
    //boil::cart.sum_real(&l2err);
    boil::cart.max_real(&linferr);
    //if(nele>0)
    //  l2err /= real(nele);
    //l2err = sqrt(l2err);
    if(n==0) {
      //l2err_first = l2err;
      linferr_first = linferr;
    }
#ifdef DEBUG
    real errrat(0);
    if(linferr_first>0) {
      errrat = linferr/linferr_first;
    }
    boil::oout<<"VOF::ev_solve: "
              <<" "<<n<<" "<<linferr<<" "<<errrat<<" ( "<<nele<<" els)"
              <<boil::endl;
#endif
    xold = x;

    /* converged? */
    //if(linferr/linferr_first<resrat) {
    if(linferr<restol) {
      converged = true;
//#ifdef DEBUG
      real errrat(0);
      if(linferr_first>0) {
        errrat = linferr/linferr_first;
      }
      if(verbose) {
        boil::oout<<"VOF::ev_solve converged after "<<n<<" steps, final abs./rel. error: "
                  <<linferr<<" "<<errrat<<"\n";
      }
//#endif
      break;
    }
  }

  if(!converged)
    boil::oout<<"VOF::ev_solve did not converge after "<<niter<<" steps, final abs./rel. error: "
                <<linferr<<" "<<linferr/linferr_first<<" !\n";

  return converged;
}
