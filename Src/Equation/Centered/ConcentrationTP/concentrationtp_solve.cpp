#include "concentrationtp.h"

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void ConcentrationTP::solve(const ResRat & fact, const char * name) {

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

  for_avijk(phi,i,j,k) {
    if(heavi->status(i,j,k)==-sig) {
      fnew[i][j][k] = phi[i][j][k];
    }
  }

  solver->solve(A, phi, fnew, MaxIter(20), name, fact);

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}
