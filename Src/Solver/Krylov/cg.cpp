#include "cg.h"

/***************************************************************************//**
*  \brief Implementation of the Conjugate Gradient (CG) solver.
*
*  \note The arguments are explained in the parent, Krylov.
*******************************************************************************/
void CG :: solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi,
                 const char * name,
                 const ResRat & res_rat, const ResTol & res_tol) {
/*----------------------------------------------+
|  templated conjugate gradient (cg) algorithm  |
+----------------------------------------------*/
  //boil::oout<<"cg:solve:begin:res_tol= "<<res_tol<<" MaxIter "<<mi<<"\n";
  //if(res_tol==boil::femto&&mi!=10){
  //  boil::oout<<"cg:point000:\n";
  //}
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
  real res = r.dot(r); 
  real res0 = res;

  // OMS(------------);
  // OPR(sqrt(res0));
  // OPR(res_tol);
  // OPR(res_rat);
 
  //if(res_tol==boil::femto&&mi!=10){
  //  boil::oout<<"cg:point100:\n";
  //}

  /* should res be scaled with A and x? */
  if(sqrt(boil::maxr(0.0,res)) < res_tol) return; // temporary meassure

  //if(res_tol==boil::femto&&mi!=10){
  //  boil::oout<<"cg:point200:\n";
  //}

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
    //if(res_tol==boil::femto&&mi!=10){
    //  boil::oout<<"cg:point300: "<<i<< "\n";
    //}
    
    /*----------+
    |  q = A p  |
    +----------*/
    q = A * p;

    /*-------------------+
    |  alfa = rho / p q  |
    +-------------------*/
    pq = p.dot(q);
    //if(res_tol==boil::femto&&mi!=10){
    //  boil::oout<<"cg:point310: "<<i<<" "<<pq<<"\n";
    //}
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
    res = r.dot(r);

    // OPR( sqrt(res) );
    //if(res_tol==boil::femto&&mi!=10){
    //  boil::oout<<"cg:point320: "<<i<<" "<<res<<"\n";
    //}

    /* should res be scaled with A and x? */
    if( sqrt(boil::maxr(0.0,res)) < res_tol ) break;

    if( sqrt(boil::maxr(0.0,res)) < sqrt(boil::maxr(0.0,res0)) * res_rat ) break; 

    /*----------------+
    |  rho_old = rho  |
    +----------------*/
    rho_old = rho;
  }
  x.exchange();

  /* for normalisation */
  q = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << sqrt(boil::maxr(0.0,res/q.dot(q))) 
                            << ", ratio = " << sqrt(boil::maxr(0.0,res/res0))
                            << ", iterations = " << i+1 
                            << boil::endl;
}
