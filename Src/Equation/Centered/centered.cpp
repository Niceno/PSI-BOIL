#include "centered.h"

/***************************************************************************//**
*  Constructor needed by coarser levels for Additive Correction method.     
*
*  \note It is recursive (i.e. calls itself).
*******************************************************************************/
Centered::Centered(const Centered * fn, const Domain * d, 
                   BndCnd & ubc, Krylov * sm) :

  Equation(d, NULL, NULL, NULL, sm), phi(*d, ubc), A(phi), 
               fnew(*d), buff(*d), res(*d), 
               fnr(fn), crsr(NULL)
{ 
  boil::oout << "Centered level " << ni()-2 << 
                            " x " << nj()-2 << 
                            " x " << nk()-2 << " created !" << boil::endl;

  if(dom->coarser() != NULL)
    crsr = new Centered( this, this->dom->coarser(), phi.bc(), solver ); 
}	

/*-----------------------------------------------------------------------------+
 '$Id: centered.cpp,v 1.17 2009/07/01 09:20:27 niceno Exp $'/
+-----------------------------------------------------------------------------*/
