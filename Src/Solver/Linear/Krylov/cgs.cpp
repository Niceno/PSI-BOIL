#include "cgs.h"

/***************************************************************************//**
*  \brief Implementation of the Conjugate Gradient Squared (CGS) solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
bool CGS :: solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi,
                  const char * name,
                  const ResRat & res_rat, const ResTol & res_tol,
                  const real scale,
                  const int stalecount,
                  const bool precform) {
/*-------------------------------------------------------+
|  templated conjugate gradient squared (cgs) algorithm  |
+-------------------------------------------------------*/

  real alfa, beta, rho, rho_old, r_v_;

  /*------------------------------------------+
  |  set the proper shape for helping arrays  | 
  +------------------------------------------*/
  p  = x.shape(); p =0.0;
  q  = x.shape(); q =0.0;
  r  = x.shape(); r =0.0;
  s  = x.shape(); s =0.0; /* used for a summ */
  u  = x.shape(); u =0.0;
  p_ = x.shape(); p_=0.0;
  q_ = x.shape(); q_=0.0;
  r_ = x.shape(); r_=0.0;
  u_ = x.shape(); u_=0.0;
  v_ = x.shape(); v_=0.0;

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

  /*------------+
  |  choose r~  |
  +------------*/
  r_ = r;

  bool converged(false);
  int i;
  for(i=0; i<mi; i++) {

    /*-------------+
    |  rho = r~ r  |
    +-------------*/
    rho = r_.dot(r);

    /*---------------------------+
    |  if rho == 0 method fails  |
    +---------------------------*/
    if( sqrt(fabs(rho)) < res_tol ) break;

    if(i == 0) {
      /*--------+
      |  u = r  |
      |  p = u  |
      +--------*/
      u = r;
      p = u;
    }
    else {
      /*----------------------------+
      |  beta = rho / rho_old       |
      |  u = r + beta q             |
      |  p = u + beta (q + beta p)  |   
      +----------------------------*/
      beta = rho / rho_old;
      u =  r + beta * q;
      p =  u + (beta*beta) * p;
      p += beta * q;
    }

    /*------------------+
    |  preconditioning  |
    +- - - - - - - - - -+
    |  solve: M p^ = p  |
    +------------------*/
    prec->solve(p_, p);

    /*------------+
    |  v^ = A p^  |
    +------------*/
    v_ = A * p_;

    /*-----------------------+
    |  alfa = rho / (r~ v^)  |
    +-----------------------*/
    r_v_ = r_.dot(v_);
    alfa = rho / r_v_;

    /*------------------+
    |  q = u - alfa v^  |
    +------------------*/
    q = u - alfa * v_;

    /*----------------------+
    |  preconditioning      |
    +- - - - - - - - - - - -+
    |  solve: M u^ = u + q  | -> should i make a new function ?
    +----------------------*/
    s = u + q;
    prec->solve(u_, s);

    /*------------------+
    |  x = x + alfa u^  |
    +------------------*/
    x += alfa * u_;

    /*------------+
    |  q^ = A u^  |
    +------------*/
    q_ = A * u_;

    /*-------------------+
    |  r  = r - alfa q^  |
    +-------------------*/
    r -= alfa * q_;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    //res = sqrt(r.dot_avg(r));
    res = sqrt(r.dot_voldiv_avg(r))/scale;

#ifdef DEBUG
    OPR(res);
#endif

    /* should res be scaled with A and u? */
    if( res < res_tol ) { converged = true; break; }

    if( res < res0 * res_rat ) { converged = true; break; }
    
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
