#include "colorcip.h"

/******************************************************************************/
void ColorCIP::totalvol() {

   /*---------+
   | method 1 |
   +---------*/
   real clrsum = 0.0;

   for_ijk(i,j,k)
     clrsum += clr[i][j][k] * dV(i,j,k);

   boil::cart.sum_real(&clrsum);

   /*---------+
   | method 2 |
   +---------*/
   /* initialize */
   for_ijk(i,j,k) {
      iflag[i][j][k]=0;
   }

   /* fluid1 */
   for_ijk(i,j,k) {
      if(clr[i][j][k]>=0.5)
      iflag[i][j][k]=1;
   }

   /* i-direction */
   for(int i=0; i<ni()-1; i++){
      for_jk(j,k){
         if((clr[i][j][k]-0.5)*(clr[i+1][j][k]-0.5)<=0.0){
            iflag[i  ][j][k]=2;
            iflag[i+1][j][k]=2;
         }
      }
   }

   /* j-direction */
   for(int j=0; j<nj()-1; j++){
      for_ik(i,k){
         if((clr[i][j][k]-0.5)*(clr[i][j+1][k]-0.5)<=0.0){
            iflag[i][j  ][k]=2;
            iflag[i][j+1][k]=2;
         }
      } 
   } 
   /* k-direction */
   for(int k=0; k<nk()-1; k++){
      for_ij(i,j){
         if((clr[i][j][k]-0.5)*(clr[i][j][k+1]-0.5)<=0.0){
            iflag[i][j][k  ]=2;
            iflag[i][j][k+1]=2;
         }
      } 
   } 

   real tvol = 0.0;

   for_ijk(i,j,k){
      if(iflag[i][j][k]==1){
         tvol += dV(i,j,k);
      } else if(iflag[i][j][k]==2) {
         tvol += clr[i][j][k] * dV(i,j,k);
      }
   }

   boil::cart.sum_real(&tvol);

   std::cout.setf(std::ios_base::scientific);
   boil::oout << "colorcip_totalvol:time,volume,clrsum= " 
              << time->current_time()
              <<" "<< tvol <<" "<< clrsum << boil::endl;
   std::cout.unsetf(std::ios_base::floatfield);

   return;
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_totalvol.cpp,v 1.3 2009/11/12 12:15:48 sato Exp $'/
+-----------------------------------------------------------------------------*/
