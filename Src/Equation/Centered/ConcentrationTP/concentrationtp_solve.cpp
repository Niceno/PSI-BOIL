#include "concentrationtp.h"

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void ConcentrationTP::solve(const ResTol & toler, const ResRat & fact,
                            const char * name) {

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
    if(heavi->status(i,j,k)==-matter_sig) {
      fnew[i][j][k] = phi[i][j][k];
    }
  }

  solver->solve(A, phi, fnew, MaxIter(20), name, fact, toler,scale*time->dti());
  phi.bnd_update();

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}
