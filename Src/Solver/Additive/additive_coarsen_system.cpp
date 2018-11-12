#include "additive.h"

/******************************************************************************/
void AC::coarsen_system(const Centered & h, Centered & H) const {
/*------------------------+ 
|  coarsen system matrix  |
+------------------------*/

  boil::timer.start("coarsening");
	
  /*---------------------------------+ 
  |  initialize central coefficient  |
  +---------------------------------*/
  for_vijk(H, iH, jH, kH) 
    H.A.c [iH][jH][kH] = 0.0;

  const int CI( (h.ni()-2) / (H.ni()-2) );
  const int CJ( (h.nj()-2) / (H.nj()-2) );
  const int CK( (h.nk()-2) / (H.nk()-2) );

  int ih[CI+1];
  int jh[CJ+1];
  int kh[CK+1];

  /*----------------------------------------+
  |  create system for coarser level of AC  |
  +----------------------------------------*/
  for_vijk(H, iH, jH, kH) {
    int iH_o = iH; // off diagonal
    int jH_o = jH; // off diagonal
    int kH_o = kH; // off diagonal

    /* create indexes */
    ih[CI] = iH * CI;
    jh[CJ] = jH * CJ;
    kh[CK] = kH * CK;
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;
    for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;
    for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

    /*--------------+
    |  c = central  |
    +--------------*/
    H.A.c[iH][jH][kH] = 0.0;
    for(int i=1; i<=CI; i++) 
      for(int j=1; j<=CJ; j++) 
        for(int k=1; k<=CK; k++) 
          H.A.c[iH][jH][kH] += h.A.c[ih[i]][jh[j]][kh[k]];

    /*------+
    |  e-w  |
    +------*/

    /* internal coefficients */
    for(int j=1; j<=CJ; j++) 
      for(int k=1; k<=CK; k++) {
        for(int i=1; i<CI; i++)
          H.A.c[iH][jH][kH] -= h.A.e[ih[i]][jh[j]][kh[k]];
        for(int i=2; i<=CI; i++) 
          H.A.c[iH][jH][kH] -= h.A.w[ih[i]][jh[j]][kh[k]];
      }

    /* external coefficients */
    H.A.e[iH_o][jH_o][kH_o] = 0.0;
    H.A.w[iH_o][jH_o][kH_o] = 0.0;
    for(int j=1; j<=CJ; j++)
      for(int k=1; k<=CK; k++) {
        H.A.e[iH_o][jH_o][kH_o] += h.A.e[ih[CI]][jh[j]][kh[k]]; 
        H.A.w[iH_o][jH_o][kH_o] += h.A.w[ih[ 1]][jh[j]][kh[k]]; 
      }

    /* restore indexes */
    jH_o = jH; 
    kH_o = kH; 
    jh[CJ] = jH * CJ;
    kh[CK] = kH * CK;
    for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;
    for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

    /*------+
    |  n-s  |
    +------*/

    /* internal coefficients */
    for(int i=1; i<=CI; i++) 
      for(int k=1; k<=CK; k++) {
        for(int j=1; j<CJ; j++)
          H.A.c[iH][jH][kH] -= h.A.n[ih[i]][jh[j]][kh[k]];
        for(int j=2; j<=CJ; j++) 
          H.A.c[iH][jH][kH] -= h.A.s[ih[i]][jh[j]][kh[k]];
      }

    /* external coefficients */
    H.A.n[iH_o][jH_o][kH_o] = 0.0;
    H.A.s[iH_o][jH_o][kH_o] = 0.0;
    for(int i=1; i<=CI; i++)
      for(int k=1; k<=CK; k++) {
        H.A.n[iH_o][jH_o][kH_o] += h.A.n[ih[i]][jh[CJ]][kh[k]]; 
        H.A.s[iH_o][jH_o][kH_o] += h.A.s[ih[i]][jh[ 1]][kh[k]]; 
      }

    /* restore indexes */
    iH_o = iH; 
    kH_o = kH; 
    ih[CI] = iH * CI;
    kh[CK] = kH * CK;
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;
    for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

    /*------+
    |  t-b  |
    +------*/

    /* internal coefficients */
    for(int i=1; i<=CI; i++) 
      for(int j=1; j<=CJ; j++) {
        for(int k=1; k<CK; k++)
          H.A.c[iH][jH][kH] -= h.A.t[ih[i]][jh[j]][kh[k]];
        for(int k=2; k<=CK; k++) 
          H.A.c[iH][jH][kH] -= h.A.b[ih[i]][jh[j]][kh[k]];
      }

    /* external coefficients */
    H.A.t[iH_o][jH_o][kH_o] = 0.0;
    H.A.b[iH_o][jH_o][kH_o] = 0.0;
    for(int i=1; i<=CI; i++)
      for(int j=1; j<=CJ; j++) {
        H.A.t[iH_o][jH_o][kH_o] += h.A.t[ih[i]][jh[j]][kh[CK]]; 
        H.A.b[iH_o][jH_o][kH_o] += h.A.b[ih[i]][jh[j]][kh[ 1]]; 
      }

  }

  /*-----------------------------------------+
  |  compute inverse of central coefficient  |
  +-----------------------------------------*/
  for_vijk(H, iH, jH, kH) 
    H.A.ci[iH][jH][kH] = 1.0 / H.A.c[iH][jH][kH];

  boil::timer.stop("coarsening");
}
