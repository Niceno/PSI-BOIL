#include "cg.h"

//#define DEBUG

/***************************************************************************//**
*  \brief Implementation of the Conjugate Gradient (CG) solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
bool CG :: solve(Matrix & A, Scalar & x, Scalar & b, const MinIter & mini,
                 const MaxIter & mi,
                 const char * name,
                 const ResRat & res_rat, const ResTol & res_tol,
                 const real scale,
                 const int stalecount,
                 const bool precform) {
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
  if(precform)
    prec->form(A, x);

  /* staleness vector */
  std::vector<real> resvect;
  if(stalecount>0) {
    resvect.resize(stalecount);
    for(auto & r : resvect)
      r = boil::unreal;
  }

  /*----------------------+
  |  compute r = b - A x  |
  +----------------------*/
  r = b - A * x;
  //real res = sqrt(r.dot_avg(r)); 
  real res = sqrt(r.dot_voldiv_avg(r))/scale;
  real res0 = res;

  if(stalecount>0) {
    resvect.back() = res;
  }

#ifdef DEBUG
  OMS(------------);
  OPR(res0);
  OPR(res_tol);
  OPR(res_rat);
#endif 

  /* should res be scaled with A and x? */
  if(res < res_tol) return true; // temporary meassure

  bool converged(false);
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
    //res = sqrt(r.dot_avg(r));
    res = sqrt(r.dot_voldiv_avg(r))/scale;

#ifdef DEBUG
    OPR(res);
    if(stalecount>0) {
      boil::oout<<i<<" ";
      for(auto & r : resvect)
        boil::oout<<" "<<r;
      boil::oout<<" "<<res<<boil::endl;
    }
#endif

    /* should res be scaled with A and x? */
    if( res < res_tol && i >= mini-1 ) { converged = true; break; }

    if( res < res0 * res_rat && i >= mini-1 ) { converged = true; break; }
    
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
        /* restore last good solution */
        //x -= alfa * p;
        break;
      } else {
        std::rotate(resvect.begin(),resvect.begin()+1,resvect.end());
        resvect.back() = res;
      }
    }

    if(i >= 30 && i%10 ==0)
      boil::oout<<"iteration= "<<i<<" residual= "<<res<<" ratio= "<<res/res0<<"\n";

    /*----------------+
    |  rho_old = rho  |
    +----------------*/
    rho_old = rho;
  } 
  x.exchange();

  /* for normalisation */
  //q = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << res
                            << ", ratio = " << res/res0
                            << ", iterations = " << i+1 
                            << boil::endl;

  return converged;
}
