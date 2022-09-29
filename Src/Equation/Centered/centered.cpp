#include "centered.h"

/***************************************************************************//**
*  Exposed constructor
*******************************************************************************/
Centered::Centered(const Domain * d,
                   const Scalar & s, const Scalar & g,
                   const Vector * v, const Times & t,
                   Matter * flu,
                   Matter * sol,
                   Linear * sm) :
  Equation(d, &t, flu, sol, sm), phi(&s), A(phi), fext(&g), u(v),
  fold(*d), fnew(*d), fbnd(*d),
  aflag(*d),
  cold(*d), cnew(*d), buff(*d),
  res(*d),
  fnr(NULL), crsr(NULL) {

  fext=phi.shape();
  fold=phi.shape();  fnew=phi.shape();  fbnd=phi.shape();
  aflag = phi.shape();
  cold=phi.shape();  cnew=phi.shape();  buff=phi.shape();
  res =phi.shape();
  set_active_flag(aflag);

  min_iter = MinIter(1);
}


/***************************************************************************//**
*  Constructor needed by coarser levels for Additive Correction method.     
*
*  The active flag is allocated but not set! That is done in AC method. 
*
*  \note It is recursive (i.e. calls itself).
*******************************************************************************/
Centered::Centered(const Centered * fn, const Domain * d, 
                   BndCnd & ubc, Linear * sm) :

  Equation(d, NULL, NULL, NULL, sm), phi(*d, ubc), A(phi), 
               fnew(*d), aflag(*d), buff(*d), res(*d), 
               fnr(fn), crsr(NULL)
{ 
  boil::oout << "Centered level "        <<
                         ni()-2*boil::BW <<
                " x " << nj()-2*boil::BW <<
                " x " << nk()-2*boil::BW << " created !" << boil::endl;

  if(dom->coarser() != NULL)
    crsr = new Centered( this, this->dom->coarser(), phi.bc(), solver ); 

  min_iter = MinIter(1);
}
