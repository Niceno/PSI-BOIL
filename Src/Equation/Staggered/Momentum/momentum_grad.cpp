#include "momentum.h"

/******************************************************************************/
void Momentum::grad(Scalar & frc) {
 
  frc.bnd_update();
  frc.exchange_all();

  /*-----------------------------------------------+
  |  add pressure gradient to the right hand side  |
  |       (initially derived from project)         |
  +-----------------------------------------------*/
  
  Comp m = Comp::u();
  for_mijk(m,i,j,k) {
#if 0
    gradp[m][i][j][k] = (frc[i-1][j][k]-frc[i][j][k]) * dSx(m,i,j,k); 
#else
    gradp[m][i][j][k] = frc[i-1][j][k]*dSx(m,Sign::neg(),i,j,k)
                      - frc[i  ][j][k]*dSx(m,Sign::pos(),i,j,k); 
#endif
  }
  
  m = Comp::v();
  for_mijk(m,i,j,k) {
#if 0
    gradp[m][i][j][k] = (frc[i][j-1][k]-frc[i][j][k]) * dSy(m,i,j,k); 
#else
    gradp[m][i][j][k] = frc[i][j-1][k]*dSy(m,Sign::neg(),i,j,k)
                      - frc[i][j  ][k]*dSy(m,Sign::pos(),i,j,k); 
#endif
  }
  
  m = Comp::w();
  for_mijk(m,i,j,k) {
#if 0
    gradp[m][i][j][k] = (frc[i][j][k-1]-frc[i][j][k]) * dSz(m,i,j,k);
#else
    gradp[m][i][j][k] = frc[i][j][k-1]*dSz(m,Sign::neg(),i,j,k)
                      - frc[i][j][k  ]*dSz(m,Sign::pos(),i,j,k); 
#endif
  }

  if(dom->ibody().nccells() > 0){
  for_m(m)
    for_mijk(m,i,j,k) 
      if( dom->ibody().off(m,i,j,k)) {
        gradp[m][i][j][k] = 0.0;
      }
  }

  return;
}
