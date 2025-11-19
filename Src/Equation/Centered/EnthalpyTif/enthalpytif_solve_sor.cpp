#include "enthalpytif.h"
using namespace std;

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void EnthalpyTIF::solve_sor(const int & it
                         , const real & omg
                         , const char * name) {

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
#if 0
  //solver->solve(A, phi, fnew, 20, name, 1e-32, fact);
  solver->solve(A, phi, fnew, MaxIter(20), name, fact);
#else
  //boil::plot->plot(*clr,phi,fnew,"clr-phi-f0", time->current_step());
  for(int mstep=1; mstep<=it; mstep++){
    for_ijk(i,j,k){
      real rhs = fnew[i][j][k]
               + A.w[i][j][k]*phi[i-1][j][k]
               + A.e[i][j][k]*phi[i+1][j][k]
               + A.s[i][j][k]*phi[i][j-1][k]
               + A.n[i][j][k]*phi[i][j+1][k]
               + A.b[i][j][k]*phi[i][j][k-1]
               + A.t[i][j][k]*phi[i][j][k+1];
      real phinew = rhs/A.c[i][j][k];
      real dphi   = phinew - phi[i][j][k];
      phi[i][j][k] += omg * dphi;
    }
    phi.exchange();
    phi.bnd_update();
  }
#endif

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}

