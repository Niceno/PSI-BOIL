#include "enthalpyfd.h"
using namespace std;

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void EnthalpyFD::solve(const ResRat & fact, const char * name) {

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

  //solver->solve(A, phi, fnew, 20, name, 1e-32, fact);
  solver->solve(A, phi, fnew, MaxIter(20), name, fact);

  phi.bnd_update();

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}

/*-----------------------------------------------------------------------------+
 '$Id: enthalpyfd_solve.cpp,v 1.4 2015/02/16 08:28:29 sato Exp $'/
+-----------------------------------------------------------------------------*/
