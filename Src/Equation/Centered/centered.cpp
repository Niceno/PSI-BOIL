#include "centered.h"

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
}
