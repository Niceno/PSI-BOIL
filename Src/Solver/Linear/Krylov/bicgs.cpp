#include "bicgs.h"

/***************************************************************************//**
*  \brief Implementation of the Bi-Conjugate Gradient Stabilized (BiCGS) 
*         solver.
*
*  \note The arguments are explained in the parent-parent, Linear.
*******************************************************************************/
void BiCGS :: solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi, 
                    const char * name,
                    const ResRat & res_rat, const ResTol & res_tol,
                    const real scale,
                    const int stalecount,
                    const bool precform) {
/*---------------------------------------------------------------+
|  templated bi-conjugate gradient stabilized (bicgs) algorithm  |
+---------------------------------------------------------------*/

  real alfa, beta, omega, rho, rho_old, r_v, ss, ts, tt;

  Range<real> small(-res_tol, res_tol);

  /*------------------------------------------+
  |  set the proper shape for helping arrays  | 
  +------------------------------------------*/
  p  = x.shape(); p =0.0;
  r  = x.shape(); r =0.0;
  s  = x.shape(); s =0.0;
  t  = x.shape(); t =0.0;
  v  = x.shape(); v =0.0;
  p_ = x.shape(); p_=0.0;
  r_ = x.shape(); r_=0.0;
  s_ = x.shape(); s_=0.0;

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
  if(res < res_tol) return; // temporary meassure

  /*------------+
  |  choose r~  |
  +------------*/
  r_ = r;

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
      |  p = r  |
      +--------*/
      p = r;
    }
    else {
      /*----------------------------------------+
      |  beta = (rho / rho_old) (alfa / omega)  |
      |  p = r + beta (beta p - omega v)        |
      +----------------------------------------*/
      beta = (rho / rho_old) * (alfa / omega);
      p =  r + beta * p;
      p -= (beta*omega) * v;
    }

    /*------------------+
    |  preconditioning  |
    +- - - - - - - - - -+
    |  solve: M p^ = p  |
    +------------------*/
    prec->solve(p_, p);

    /*-----------+
    |  v = A p^  |
    +-----------*/
    v = A * p_;

    /*--------------------+
    |  alfa = rho / r~ v  |
    +--------------------*/
    r_v = r_.dot(v);
    alfa = rho / r_v;

    /*-----------------+
    |  s = r - alfa v  |
    +-----------------*/
    s = r - alfa * v;

    /*------------------------------------+
    |  chack norm of s; if small enough,  |
    |    set x = x + alfa p^ and stop     |
    +------------------------------------*/
    ss = s.dot(s);
    if( sqrt(ss) < res_tol ) {
      x += alfa * p_;
      break;
    }

    /*------------------+
    |  preconditioning  |
    +- - - - - - - - - -+
    |  solve: M s^ = s  |
    +------------------*/
    prec->solve(s_, s);

    /*-------------+ 
    |  t = A * s^  |
    +-------------*/
    t = A * s_;

    /*----------------------+
    |  omega = (t s)/(t t)  |
    +----------------------*/
    ts = t.dot(s);
    tt = t.dot(t);
    omega = ts / tt;

    /*-----------------------------+
    |  x = x + alfa p^ + omega s^  |
    +-----------------------------*/
    x += alfa  * p_;
    x += omega * s_;

    /*-------------------+
    |  r  = s - omega t  |
    +-------------------*/
    r = s - omega * t;

    /*--------------------+
    |  exit if converged  |
    +--------------------*/
    //res = sqrt(r.dot_avg(r));
    res = sqrt(r.dot_voldiv_avg(r))/scale;

#ifdef DEBUG
    OPR(res);
#endif

    /* should res be scaled with A and u? */
    if( res < res_tol ) break;

    if( res < res0 * res_rat ) break; 

    /*---------------------------------------------------+
    |  for continuation it is necessary that omega != 0  |
    +---------------------------------------------------*/
    if( small.contains(omega) ) break;

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
  //p = A * x;

  if(name!=NULL) boil::oout << name 
                            << ", residual = " << res
                            << ", ratio = " << res/res0
                            << ", iterations = " << i+1 
                            << boil::endl;
}
