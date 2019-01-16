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

  const int CI( (h.ei()-h.si()+1) / (H.ei()-H.si()+1) );
  const int CJ( (h.ej()-h.sj()+1) / (H.ej()-H.sj()+1) );
  const int CK( (h.ek()-h.sk()+1) / (H.ek()-H.sk()+1) );

  int ih[CI+2], jh[CJ+2], kh[CK+2];

  for(int iH=H.si(); iH<=H.ei(); iH++) {

    /* create indexes */
    ih[CI] = iH*CI - (H.si() - 1) * (CI - 1);
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;

    for(int jH=H.sj(); jH<=H.ej(); jH++) {

      /* create indexes */
      jh[CJ] = jH*CJ - (H.sj() - 1) * (CJ - 1);
      for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;

      for(int kH=H.sk(); kH<=H.ek(); kH++) {

        /* create indexes */
        kh[CK] = kH*CK - (H.sk() - 1) * (CK - 1);
        for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++)
            for(int k=1; k<=CK; k++)
              h.phi[ih[i]][jh[j]][kh[k]] += H.phi[iH][jH][kH];
      }
    }
  }
}
