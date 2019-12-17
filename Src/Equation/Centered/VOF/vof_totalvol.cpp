#include "vof.h"

/******************************************************************************/
void VOF::totalvol() {

   /*---------+
   | method 1 |
   +---------*/
   real phisum = 0.0;
   real vapsum = 0.0;

   for_ijk(i,j,k){
     if(dom->ibody().off(i,j,k)) continue;
     phisum += phi[i][j][k] * dV(i,j,k);
     vapsum += (1.0-phi[i][j][k]) * dV(i,j,k);
   }

   boil::cart.sum_real(&phisum);
   boil::cart.sum_real(&vapsum);


   std::cout.setf(std::ios_base::scientific);
   boil::oout << "totalvol:time,volume,phisum= " 
              << time->current_time()
              <<" "<< phisum << " "<<vapsum<< boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);

   return;
}
