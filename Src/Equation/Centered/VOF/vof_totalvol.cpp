#include "vof.h"

/******************************************************************************/
void VOF::totalvol() {

   /*---------+
   | method 1 |
   +---------*/
   real phisum = 0.0;
   real phisum_inv = 0.0;

   for_ijk(i,j,k){
     if(dom->ibody().off(i,j,k)) continue;
     phisum += phi[i][j][k] * dV(i,j,k);
     phisum_inv += (1.0-phi[i][j][k]) * dV(i,j,k);
   }

   boil::cart.sum_real(&phisum);
   boil::cart.sum_real(&phisum_inv);

   std::cout.setf(std::ios_base::scientific);
   std::cout<< std::setprecision(16);
   boil::oout << "totalvol:time,volume,phisum= " 
              << time->current_time()
              <<" "<< phisum << " "<< phisum_inv << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);
   std::cout<< std::setprecision(6);

   total_vol0 = phisum;
   total_vol1 = phisum_inv;

   return;
}
