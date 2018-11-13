#include "enthalpytif.h"
using namespace std;

/***************************************************************************//**
*  Solves linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*******************************************************************************/
void EnthalpyTIF::deltat(Scalar & deltaT, const Scalar & heaviside,
                         const ResRat & fact, const char * name, const real factor) {

  assert(solver);

  tint_field(heaviside,factor,true);

  /* set the name and start timing */
  std::string msg;
  if( name ) 
   {
    msg = name;
    msg += " solver";
    
    boil::timer.start(msg.c_str());
   }

  /* solve */
  ftifold = ftif;
  update_ftif();
  fdelta = (ftif - ftifold);
  //fdelta *= 0.01;
 
#if 0
  for_aijk(i,j,k) {
    if(j==1&&k==1) boil::oout<<" "<<i<<" "<<phi[i][j][k]-tsat<<" "<<tif[i][j][k]-tsat<<" "<<deltaT[i][j][k]<<" "<<ftif[i][j][k]<<" "<<ftifold[i][j][k]<<boil::endl;
  }
#endif

  solver->solve(A, deltaT, fdelta, MaxIter(20), name, fact);

  deltaT.bnd_update();

  //boil::oout << "Here! "<< tif[50][1][1]<<" "<<phi[50][1][1]<<" "<<phi[51][1][1]<<boil::endl;
#if 0
  for_aijk(i,j,k) {
    if(j==1&&k==1) boil::oout<<" "<<i<<" "<<phi[i][j][k]-tsat<<" "<<tif[i][j][k]-tsat<<" "<<deltaT[i][j][k]<<" "<<fdelta[i][j][k]<<boil::endl;
  }
#endif

  /* stop the timing */
  if( name ) 
   {boil::timer.stop(msg.c_str());}
}

