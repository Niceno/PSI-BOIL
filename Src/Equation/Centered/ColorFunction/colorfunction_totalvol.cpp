#include "colorfunction.h"
#include <iomanip>

/******************************************************************************/
void ColorFunction::totalvol() {
/***************************************************************************//**
*  \brief Calculate volume by two methods.
*         Method1: Sum(phi*dV)
*         Method2: Sum(dV)  in the region phi > phisurf
*******************************************************************************/

   clrsum1=0.0;
   clrsum2=0.0;

   if(dom->ibody().nccells()==0){
     for_ijk(i,j,k){
       clrsum1 += phi[i][j][k] * dV(i,j,k);
       clrsum2 += (1.0-phi[i][j][k]) * dV(i,j,k);
     }
   } else {
     for_ijk(i,j,k){
       clrsum1 += phi[i][j][k]       * dV(i,j,k) * dom->ibody().fV(i,j,k);
       clrsum2 += (1.0-phi[i][j][k]) * dV(i,j,k) * dom->ibody().fV(i,j,k);
     }
   }

   boil::cart.sum_real(&clrsum1);
   boil::cart.sum_real(&clrsum2);

#if 0
   std::cout.setf(std::ios_base::scientific);
   std::cout<< std::setprecision(16);
   boil::oout << "cipcsl2_totalvol:time,clrsum1,clrsum2= " 
              << time->current_time()
              <<" "<< clrsum1 <<" "<< clrsum2 << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);
   std::cout<< std::setprecision(6);
#endif

   return;
}
