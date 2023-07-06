#include "enthalpyfd.h"

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void EnthalpyFD::solve(const ResTol & toler, const ResRat & fact, 
                       const MaxIter & maxiter, const char * name) {

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

  //solver->solve(A, phi, fnew, MinIter(1), MaxIter(1000), 
  //solver->solve(A, phi, fnew, min_iter, MaxIter(1000), 
  solver->solve(A, phi, fnew, min_iter, maxiter, 
                name, fact, toler,scale*time->dti());

  phi.bnd_update();

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}
