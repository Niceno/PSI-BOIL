#include "centered.h"

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void Centered::solve(const ResTol & toler, const ResRat & fact, 
                     const char * name, int itmax) {

  assert(solver);

  /* set the name and start timing */
  std::string msg;
  if( name ) 
   {
    msg = name;
    msg += " solver";
    
    boil::timer.start(msg.c_str());
   }

  /* solve */
  update_rhs();

  //if (itmax==NULL) { int itmp = 20; itmax = &itmp; }
  //const int it = *itmax;
  if (itmax==NULL) { itmax = 20; }
  const int it = itmax;
  solver->solve(A, phi, fnew, min_iter, MaxIter(it), name, fact, toler,scale*time->dti());
  phi.bnd_update();

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}
