#include "centered.h"

/******************************************************************************/
void Centered::coarsen() 
{
  OMS(coarsenig);

  if(dom->coarser() != NULL && crsr == NULL)
    crsr = new Centered( this, this->dom->coarser(), phi.bc(), solver ); 

//if(this->coarser() != NULL)
//  crsr->coarsen_system();
}

/*-----------------------------------------------------------------------------+
 '$Id: centered_coarsen.cpp,v 1.7 2009/07/01 09:20:27 niceno Exp $'/
+-----------------------------------------------------------------------------*/
