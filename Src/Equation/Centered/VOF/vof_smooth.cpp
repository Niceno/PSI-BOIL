#include "vof.h"

/******************************************************************************/
void VOF::smooth(){
/******************************************************************************/
/* Smooth volume fraction to minimize parasitic current. (Benchmark numerical */
/* simulations of segmented two-phase flows in microchannels using the Volume */
/* of Fluid method)                                                           */
/******************************************************************************/
  for(int layer=1; layer<=2; layer++){
    for_ijk(i,j,k){
      stmp2[i][j][k] = ((stmp[i][j][k] + stmp[i-1][j][k])/2.0 * dSx(i,j,k) +
                        (stmp[i][j][k] + stmp[i+1][j][k])/2.0 * dSx(i,j,k) +
                        (stmp[i][j][k] + stmp[i][j-1][k])/2.0 * dSy(i,j,k) +
                        (stmp[i][j][k] + stmp[i][j+1][k])/2.0 * dSy(i,j,k) +
                        (stmp[i][j][k] + stmp[i][j][k-1])/2.0 * dSz(i,j,k) +
                        (stmp[i][j][k] + stmp[i][j][k+1])/2.0 * dSz(i,j,k))/
                        (2 * (dSx(i,j,k)+dSy(i,j,k)+dSz(i,j,k)));
    } 
    stmp2.exchange();
    for_aijk(i,j,k){
      stmp[i][j][k] = stmp2[i][j][k];
    }
    stmp.exchange();
  }
  return;                
}  
