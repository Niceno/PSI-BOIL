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
