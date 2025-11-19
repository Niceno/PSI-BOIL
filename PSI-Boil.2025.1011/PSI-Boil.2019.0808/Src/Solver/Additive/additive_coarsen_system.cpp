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

  const int CI( (h.ei()-h.si()+1) / (H.ei()-H.si()+1) );
  const int CJ( (h.ej()-h.sj()+1) / (H.ej()-H.sj()+1) );
  const int CK( (h.ek()-h.sk()+1) / (H.ek()-H.sk()+1) );

  int ih[CI+1];
  int jh[CJ+1];
  int kh[CK+1];

  /*----------------------------------------+
  |  create system for coarser level of AC  |
  +----------------------------------------*/
  for_vijk(H, iH, jH, kH) {

    /* create indexes */
    ih[CI] = iH*CI - (H.si() - 1) * (CI - 1);
    jh[CJ] = jH*CJ - (H.sj() - 1) * (CJ - 1);
    kh[CK] = kH*CK - (H.sk() - 1) * (CK - 1);
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;
    for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;
    for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

    /* active cell flagging */
    int central_flag = 0;

    /* if cell is inactive, central flag must be 1
       to avoid singularity */
    if(!H.aflag[iH][jH][kH]) {
      central_flag = 1;
    }

    /* w,e,s,n,b,t elements are zero for inactive cells */

    /*--------------+
    |  c = central  |
    +--------------*/
    H.A.c[iH][jH][kH] = 0.0;
    for(int i=1; i<=CI; i++) 
      for(int j=1; j<=CJ; j++) 
        for(int k=1; k<=CK; k++) {
          H.A.c[iH][jH][kH] += h.A.c[ih[i]][jh[j]][kh[k]] 
#if 0 /* discard inactive cells from blending (sorry for the syntax) */ 
                             * std::max(central_flag,
                                        h.aflag[ih[i]][jh[j]][kh[k]])
#endif
                             ;
        }

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
    H.A.e[iH][jH][kH] = 0.0;
    H.A.w[iH][jH][kH] = 0.0;
    for(int j=1; j<=CJ; j++)
      for(int k=1; k<=CK; k++) {
        H.A.e[iH][jH][kH] += h.A.e[ih[CI]][jh[j]][kh[k]]; 
        H.A.w[iH][jH][kH] += h.A.w[ih[ 1]][jh[j]][kh[k]]; 
      }

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
    H.A.n[iH][jH][kH] = 0.0;
    H.A.s[iH][jH][kH] = 0.0;
    for(int i=1; i<=CI; i++)
      for(int k=1; k<=CK; k++) {
        H.A.n[iH][jH][kH] += h.A.n[ih[i]][jh[CJ]][kh[k]]; 
        H.A.s[iH][jH][kH] += h.A.s[ih[i]][jh[ 1]][kh[k]]; 
      }

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
    H.A.t[iH][jH][kH] = 0.0;
    H.A.b[iH][jH][kH] = 0.0;
    for(int i=1; i<=CI; i++)
      for(int j=1; j<=CJ; j++) {
        H.A.t[iH][jH][kH] += h.A.t[ih[i]][jh[j]][kh[CK]]; 
        H.A.b[iH][jH][kH] += h.A.b[ih[i]][jh[j]][kh[ 1]]; 
      }

  }

  /*-----------------------------------------+
  |  compute inverse of central coefficient  |
  +-----------------------------------------*/
  for_vijk(H, iH, jH, kH) {
    H.A.ci[iH][jH][kH] = 1.0 / H.A.c[iH][jH][kH];
  }

  boil::timer.stop("coarsening");
}
