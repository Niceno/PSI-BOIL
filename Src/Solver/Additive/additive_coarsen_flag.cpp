#include "additive.h"

/******************************************************************************/
void AC::coarsen_flag(const Centered & h, Centered & H) const {
/*----------------------+ 
|  coarsen active flag  |
+----------------------*/

  const int CI( (h.ei()-h.si()+1) / (H.ei()-H.si()+1) );
  const int CJ( (h.ej()-h.sj()+1) / (H.ej()-H.sj()+1) );
  const int CK( (h.ek()-h.sk()+1) / (H.ek()-H.sk()+1) );

  int ih[CI+1];
  int jh[CJ+1];
  int kh[CK+1];

  /*-----------------------------------------+
  |  calculate flag for coarser level of AC  |
  +-----------------------------------------*/
  for_vijk(H, iH, jH, kH) {

    /* create indexes */
    ih[CI] = iH*CI - (H.si() - 1) * (CI - 1);
    jh[CJ] = jH*CJ - (H.sj() - 1) * (CJ - 1);
    kh[CK] = kH*CK - (H.sk() - 1) * (CK - 1);
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;
    for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;
    for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

    /* simply use 1 if any flag on fine level is 1 */
    H.aflag[iH][jH][kH] = 0;
    for(int i=1; i<=CI; i++) 
      for(int j=1; j<=CJ; j++) 
        for(int k=1; k<=CK; k++) {
          /* compound boolean-or and assignment operator */
          H.aflag[iH][jH][kH] |= h.aflag[ih[i]][jh[j]][kh[k]];
        }
  }

  return;
}
