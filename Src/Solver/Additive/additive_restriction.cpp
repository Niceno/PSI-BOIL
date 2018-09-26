#include "additive.h"

/***************************************************************************//**
*  \param h - Centered variable on fine grid,
*  \param H - Centered variable on coarse grid.
*
*  \note The algorithm is flexible and resolution ratio between fine and 
*        coarse grid does not have to be two. It can be any integer value. 
*        That is why grid resolution do not need number of cells to be a
*        power of two.
*******************************************************************************/
void AC::restriction(const Centered & h, Centered & H) const {
/*-----------------------------------------------------------------+
|                                                                  |
|                   +-------------------------+                    |
|                   |  R E S T R I C T I O N  |                    |
|                   +-------------------------+                    |
|                                                                  |
|                                                                  |
|   from fine to coarse:                                           |
|   ~~~~~~~~~~~~~~~~~~~~                                           |
|                                                                  |
|   h:    | . |-O-|-O-|-O-|-O-|-O-|-O-|-O-|-O-| . |    h.ni = 10   |
|           0   1   2   3   4   5   6   7   8   9                  |
|                                                                  |
|                                                                  |
|   H:  |  .  |---O---|---O---|---O---|---O---|  .  |  H.ni =  6   |
|          0      1       2       3       4      5                 |
|                                                                  |
|                                                                  |
+-----------------------------------------------------------------*/

  const int CI( (h.ni()-2) / (H.ni()-2) );
  const int CJ( (h.nj()-2) / (H.nj()-2) );
  const int CK( (h.nk()-2) / (H.nk()-2) );
  int ih_[CI+2], jh_[CJ+2], kh_[CK+2];
  int ih [CI+2], jh [CJ+2], kh [CK+2];

  real D[CI+1][CJ+1][CK+1]; // [0] not used
	
  h.phi.exchange();

  for(int iH=1; iH < H.ni()-1; iH++) {
        /* create indexes */
        ih_[CI]   = iH * CI;
        ih_[CI+1] = ih_[CI] + 1;
        for(int i=1; i<=CI; i++) ih_[CI-i] = ih_[CI]-i;

    for(int jH=1; jH < H.nj()-1; jH++) {
          /* create indexes */
          jh_[CJ]   = jH * CJ;
          jh_[CJ+1] = jh_[CJ] + 1;
          for(int j=1; j<=CJ; j++) jh_[CJ-j] = jh_[CJ]-j;
    
      for(int kH=1; kH < H.nk()-1; kH++) {
            /* create indexes */
            kh_[CK]   = kH * CK;
            kh_[CK+1] = kh_[CK] + 1;
            for(int k=1; k<=CK; k++) kh_[CK-k] = kh_[CK]-k;

        /*-----------------------------------------------------------------+ 
        |  compute new diagonal terms and start computing right hand side  |
        +-----------------------------------------------------------------*/ 
        H.fnew[iH][jH][kH] = 0.0;
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++)
            for(int k=1; k<=CK; k++) {
              D[i][j][k]          = h.A.c [ih_[i]][jh_[j]][kh_[k]];
              H.fnew[iH][jH][kH] += h.fnew[ih_[i]][jh_[j]][kh_[k]];
            }

        /* update indexes */
        for(int i=1; i<=CI; i++) ih[i] = ih_[i];
        for(int j=1; j<=CJ; j++) jh[j] = jh_[j];
        for(int k=1; k<=CK; k++) kh[k] = kh_[k];

        /*------+
        |  e-w  |
        +------*/
        for(int j=1; j<=CJ; j++)
          for(int k=1; k<=CK; k++) {
            H.fnew[iH][jH][kH] += h.A.e[ih_[CI  ]][jh [j]][kh [k]] 
                                * h.phi[ih_[CI+1]][jh_[j]][kh_[k]];
            H.fnew[iH][jH][kH] += h.A.w[ih_[1   ]][jh [j]][kh [k]] 
                                * h.phi[ih_[0   ]][jh_[j]][kh_[k]];

            for(int i=1; i< CI; i++)
              D[i][j][k] -= (h.A.e[ih_[i]][jh[j]][kh[k]]);

            for(int i=CI; i>1; i--)
              D[i][j][k] -= (h.A.w[ih_[i]][jh[j]][kh[k]]);
          }

        /* update indexes */
        for(int j=1; j<=CJ; j++) jh[j] = jh_[j];
        for(int k=1; k<=CK; k++) kh[k] = kh_[k];

        /*------+
        |  n-s  |
        +------*/
        for(int i=1; i<=CI; i++)
          for(int k=1; k<=CK; k++) {
            H.fnew[iH][jH][kH] += h.A.n[ih [i]][jh_[CJ  ]][kh [k]] 
                                * h.phi[ih_[i]][jh_[CJ+1]][kh_[k]];
            H.fnew[iH][jH][kH] += h.A.s[ih [i]][jh_[1   ]][kh [k]] 
                                * h.phi[ih_[i]][jh_[0   ]][kh_[k]];

            for(int j=1; j< CJ; j++)
              D[i][j][k] -= (h.A.n[ih[i]][jh_[j]][kh[k]]);

            for(int j=CJ; j>1; j--)
              D[i][j][k] -= (h.A.s[ih[i]][jh_[j]][kh[k]]);
          }

        /* update indexes */
        for(int i=1; i<=CI; i++) ih[i] = ih_[i];
        for(int k=1; k<=CK; k++) kh[k] = kh_[k];

        /*------+
        |  t-b  |
        +------*/
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++) {
            H.fnew[iH][jH][kH] += h.A.t[ih [i]][jh [j]][kh_[CK  ]]
                                * h.phi[ih_[i]][jh_[j]][kh_[CK+1]];
            H.fnew[iH][jH][kH] += h.A.b[ih [i]][jh [j]][kh_[1   ]]
                                * h.phi[ih_[i]][jh_[j]][kh_[0   ]];

            for(int k=1; k< CK; k++)
              D[i][j][k] -= h.A.t[ih[i]][jh[j]][kh_[k]];

            for(int k=CK; k> 1; k--)
              D[i][j][k] -= h.A.b[ih[i]][jh[j]][kh_[k]];
          }

        /*---------------------------+ 
        |  finalize right hand side  |
        +---------------------------*/ 
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++) 
            for(int k=1; k<=CK; k++) 
              H.fnew[iH][jH][kH] -= D[i][j][k] * h.phi[ih_[i]][jh_[j]][kh_[k]];
      }	
    }  
  }	  

  /*----------------------------------------------------+ 
  |  r.h.s. normalization.                              |
  + - - - - - - - - - - - - - - - - - - - - - - - - - - +
  |  maybe this is quite normal. a.c. should solve      |
  |  coarse levels directly. if i do not do that, i     |
  |  get distorted r.h.s. on even coarser levels        |
  +----------------------------------------------------*/
  real sum = 0.0;
  for(int i=1; i<=H.ni()-2; i++)
    for(int j=1; j<=H.nj()-2; j++)
      for(int k=1; k<=H.nk()-2; k++) 
        sum += H.fnew[i][j][k];

  boil::cart.sum_real(& sum);

  int ncell = (H.ni()-2)*(H.nj()-2)*(H.nk()-2);
  boil::cart.sum_int( & ncell );

  sum /= (real)ncell;

  for(int i=1; i<=H.ni()-2; i++)
    for(int j=1; j<=H.nj()-2; j++)
      for(int k=1; k<=H.nk()-2; k++) 
        H.fnew[i][j][k] -= sum;
}

/*-----------------------------------------------------------------------------+
 '$Id: additive_restriction.cpp,v 1.11 2009/11/03 12:25:38 niceno Exp $'/
+-----------------------------------------------------------------------------*/
