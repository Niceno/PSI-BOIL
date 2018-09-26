#include "centered.h"

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void Centered::solve(const ResRat & fact, const char * name) {

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
  solver->solve(A, phi, fnew, MaxIter(20), name, fact);

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}

/*-----------------------------------------------------------------------------+
 '$Id: centered_solve.cpp,v 1.2 2011/05/25 12:29:01 niceno Exp $'/
+-----------------------------------------------------------------------------*/
