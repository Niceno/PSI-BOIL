#include "levelset.h"

/******************************************************************************/
void LevelSet::totalvol(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate sum of sca
*******************************************************************************/
   real scasum1 = 0.0;
   real scasum2 = 0.0;

   for_ijk(i,j,k){
     scasum1 += sca[i][j][k] * dV(i,j,k);
     scasum2 += (1.0-sca[i][j][k]) * dV(i,j,k);
  }
   boil::cart.sum_real(&scasum1);
   boil::cart.sum_real(&scasum2);

   std::cout.setf(std::ios_base::scientific);
   boil::oout << "levelset_totalvol:time,sum1,sum2= " 
              << time->current_time()
              <<" "<< scasum1 <<" "<< scasum2 << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);

   return;
}
