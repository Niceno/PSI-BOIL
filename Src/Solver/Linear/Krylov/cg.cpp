#include "cg.h"

//#define DEBUG

/***************************************************************************//**
*  \brief Implementation of the Conjugate Gradient (CG) solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
void CG :: solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi,
                 const char * name,
                 const ResRat & res_rat, const ResTol & res_tol) {
/*----------------------------------------------+
|  templated conjugate gradient (cg) algorithm  |
+----------------------------------------------*/
  real alfa, beta, rho, rho_old, pq;

  /*------------------------------------------+
  |  set the proper shape for helping arrays  | 
  +------------------------------------------*/
  p = x.shape(); p=0.0;
  q = x.shape(); q=0.0;
  r = x.shape(); r=0.0;
  z = x.shape(); z=0.0;

  /*--------------------------------+
  |  form preconditioning matrix M  |
  +--------------------------------*/
  prec->form(A, x);

  /*----------------------+
  |  compute r = b - A x  |
  +----------------------*/
  r = b - A * x;
  //real res = sqrt(r.dot(r)); 
  real res = sqrt(r.dot_voldiv_avg(r));
  real res0 = res;

#ifdef DEBUG
  OMS(------------);
  OPR(res0);
  OPR(res_tol);
  OPR(res_rat);
#endif

  /* should res be scaled with A and x? */
  if(res < res_tol) return; // temporary meassure

  int i;
  for(i=0; i<mi; i++) {

    /*-----------------+
    |  preconditioning |
    +- - - - - - - - - +
    |  solve: M z = r  |
    +------------------*/
    prec->solve(z, r);
  
    /*------------+
    |  rho = r z  |
    +------------*/
    rho = r.dot(z);
	
    if(i == 0) {
      /*--------+    
      |  p = z  |
      +--------*/    
      p = z;
      
    } else {
      /*-----------------------+
      |  beta = rho / rho_old  |
      |  p = z + beta * p      |
      +-----------------------*/
      beta = rho / rho_old;
      p = z + beta * p;
    }
    //p.exchange(); //unnecessary because exchange will be carried out when
    // A * p is called. It is implemented in scalar_operators.cpp //
    
    /*----------+
    |  q = A p  |
    +----------*/
    q = A * p;

    /*-------------------+
    |  alfa = rho / p q  |
    +-------------------*/
    pq = p.dot(q);
    alfa = rho / pq;
    
    /*-----------------+ 
    |  x = x + alfa p  |
    |  r = r - alfa q  |
    +-----------------*/ 
    x += alfa * p;
    r -= alfa * q;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    //res = sqrt(r.dot(r));
    res = sqrt(r.dot_voldiv_avg(r));

#ifdef DEBUG
    OPR(res);
#endif

    /* should res be scaled with A and x? */
    if( res < res_tol ) break;

    if( res < res0 * res_rat ) break; 

    /*----------------+
    |  rho_old = rho  |
    +----------------*/
    rho_old = rho;
  }
  x.exchange();

  /* for normalisation */
  q = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << res/sqrt(q.dot_voldiv_avg(q))
                            << ", ratio = " << res/res0
                            << ", iterations = " << i+1 
                            << boil::endl;
}
