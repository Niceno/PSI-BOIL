#include "tif.h"

/***************************************************************************//**
 *  Calculates interface temperature field
 *  Tn+1 = Tint*factor + Tn*(1-factor)
*******************************************************************************/
void TIF::extend_tint() { /* crude code */
  for_vijk(tif,i,j,k) {
    if(tif.domain()->ibody().off(i,j,k)) continue;
    if(!Interface(i,j,k) && Vicinity(i,j,k)) {
      int ctr(0);
      real tintf(0.0);

      if(tif.domain()->ibody().on(i+1,j,k)&&Interface(i+1,j,k)) {
        ctr++;
        tintf += tif[i+1][j][k];
      } 
      if(tif.domain()->ibody().on(i-1,j,k)&&Interface(i-1,j,k)) {
        ctr++;
        tintf += tif[i-1][j][k];
      } 
      if(tif.domain()->ibody().on(i,j+1,k)&&Interface(i,j+1,k)) {
        ctr++;
        tintf += tif[i][j+1][k];
      } 
      if(tif.domain()->ibody().on(i,j-1,k)&&Interface(i,j-1,k)) {
        ctr++;
        tintf += tif[i][j-1][k];
      } 
      if(tif.domain()->ibody().on(i,j,k+1)&&Interface(i,j,k+1)) {
        ctr++;
        tintf += tif[i][j][k+1];
      } 
      if(tif.domain()->ibody().on(i,j,k-1)&&Interface(i,j,k-1)) {
        ctr++;
        tintf += tif[i][j][k-1];
      }  

      if(ctr > 0) {
        tintf /= real(ctr);
        tif[i][j][k] = tintf;
      }
    }
  }

  tif.bnd_update();
  tif.exchange_all();

  return;
}

