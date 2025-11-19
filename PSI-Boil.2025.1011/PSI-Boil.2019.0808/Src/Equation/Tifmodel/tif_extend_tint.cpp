#include "tif.h"

/***************************************************************************//**
 *  Calculates interface temperature field
 *  Tn+1 = Tint*factor + Tn*(1-factor)
*******************************************************************************/
void TIF::extend_tint() { /* crude code */

  for_avijk(tif,i,j,k) {
    tempflag[i][j][k] = abs(iflag[i][j][k])<2 ? 1 : 0;
  }

  stmp = tif;
  tempflag2 = tempflag;
  
#if 1
  for(int iloop=1; iloop<3; iloop++) { 
    for_vijk(tif,i,j,k) {
      if(tempflag[i][j][k]==0) {
  
        int iflagw = tempflag[i-1][j][k];
        int iflage = tempflag[i+1][j][k];
        int iflags = tempflag[i][j-1][k];
        int iflagn = tempflag[i][j+1][k];
        int iflagb = tempflag[i][j][k-1];
        int iflagt = tempflag[i][j][k+1];

        int inb = iflagw+iflage+iflags+iflagn+iflagb+iflagt;

        if(inb>0) {
          stmp[i][j][k] = real(  iflagw*tif[i-1][j][k] + iflage*tif[i+1][j][k]
                               + iflags*tif[i][j-1][k] + iflagn*tif[i][j+1][k]
                               + iflagb*tif[i][j][k-1] + iflagt*tif[i][j][k+1]
                              )
                        / real(inb);

          tempflag2[i][j][k] = 2; /* tempflag=2 for extrapolated */
        }
      }
    }
    stmp.bnd_update_symmetry(); /* copy on symmetry plane */
    tempflag2.bnd_update_symmetry(); /* copy on symmetry plane */
    stmp.exchange();
    tempflag2.exchange();
    tif = stmp;
    tempflag = tempflag2;
  }
#endif

  tif.bnd_update();
  tif.exchange_all();

  return;
}

