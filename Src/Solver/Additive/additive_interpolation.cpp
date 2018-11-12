#include "additive.h"

/***************************************************************************//**
*  \param H - Centered variable on coarse grid,
*  \param h - Centered variable on fine grid.
*
*  \note The algorithm is flexible and resolution ratio between fine and 
*        coarse grid does not have to be two. It can be any integer value. 
*        That is why grid resolution do not need number of cells to be a
*        power of two.
*******************************************************************************/
void AC::interpolation(const Centered & H, Centered & h) const { 
/*-----------------------------------------------------------------+
|                                                                  |
|                 +-----------------------------+                  |
|                 |  I N T E R P O L A T I O N  |                  |
|                 +-----------------------------+                  |
|                                                                  |
|                                                                  |
|   from coarse to fine:  (from bottom to top)                     |
|   ~~~~~~~~~~~~~~~~~~~~                                           |
|                                                                  |
|                                                                  |
|   h:    | . |-O-|-O-|-O-|-O-|-O-|-O-|-O-|-O-| . |    h.ni = 10   |
|           0   1   2   3   4   5   6   7   8   9                  |
|                                                                  |
|                                                                  |
|   H:  |  .  |---O---|---O---|---O---|---O---|  .  |  H.ni =  6   |
|          0      1       2       3       4      5                 |
|                                                                  |
+-----------------------------------------------------------------*/

  const int CI( (h.ni()-2) / (H.ni()-2) );
  const int CJ( (h.nj()-2) / (H.nj()-2) );
  const int CK( (h.nk()-2) / (H.nk()-2) );

  int ih_[CI+2], jh_[CJ+2], kh_[CK+2];

  for(int iH=1; iH < H.ni()-1; iH++) {
        /* create indexes */
        ih_[CI]   = iH * CI;
        for(int i=1; i<=CI; i++) ih_[CI-i] = ih_[CI]-i;

    for(int jH=1; jH < H.nj()-1; jH++) {
          /* create indexes */
          jh_[CJ]   = jH * CJ;
          for(int j=1; j<=CJ; j++) jh_[CJ-j] = jh_[CJ]-j;

      for(int kH=1; kH < H.nk()-1; kH++) {
            /* create indexes */
            kh_[CK]   = kH * CK;
            for(int k=1; k<=CK; k++) kh_[CK-k] = kh_[CK]-k;

        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++)
            for(int k=1; k<=CK; k++)
              h.phi[ih_[i]][jh_[j]][kh_[k]] += H.phi[iH][jH][kH];
      }  
    }  
  }	  
}
