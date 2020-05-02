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

  const int CI( (h.ei()-h.si()+1) / (H.ei()-H.si()+1) );
  const int CJ( (h.ej()-h.sj()+1) / (H.ej()-H.sj()+1) );
  const int CK( (h.ek()-h.sk()+1) / (H.ek()-H.sk()+1) );

  int ih[CI+2], jh[CJ+2], kh[CK+2];

  real D[CI+1][CJ+1][CK+1]; // [0] not used
	
  h.phi.exchange();

  for(int iH=H.si(); iH<=H.ei(); iH++) {

    /* create indexes */
    ih[CI]   = iH*CI - (H.si() - 1) * (CI - 1);
    ih[CI+1] = ih[CI] + 1;
    for(int i=1; i<=CI; i++) ih[CI-i] = ih[CI]-i;

    for(int jH=H.sj(); jH<=H.ej(); jH++) {

      /* create indexes */
      jh[CJ]   = jH*CJ - (H.sj() - 1) * (CJ - 1);
      jh[CJ+1] = jh[CJ] + 1;
      for(int j=1; j<=CJ; j++) jh[CJ-j] = jh[CJ]-j;

      for(int kH=H.sk(); kH<=H.ek(); kH++) {

        /* create indexes */
        kh[CK]   = kH*CK - (H.sk() - 1) * (CK - 1);
        kh[CK+1] = kh[CK] + 1;
        for(int k=1; k<=CK; k++) kh[CK-k] = kh[CK]-k;

        /*-----------------------------------------------------------------+ 
        |  compute new diagonal terms and start computing right hand side  |
        +-----------------------------------------------------------------*/ 
        H.fnew[iH][jH][kH] = 0.0;
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++)
            for(int k=1; k<=CK; k++) {
              D[i][j][k]          = h.A.c [ih[i]][jh[j]][kh[k]];
              H.fnew[iH][jH][kH] += h.fnew[ih[i]][jh[j]][kh[k]];
            }

        /*------+
        |  e-w  |
        +------*/
        for(int j=1; j<=CJ; j++)
          for(int k=1; k<=CK; k++) {
            H.fnew[iH][jH][kH] += h.A.e[ih[CI  ]][jh[j]][kh[k]] 
                                * h.phi[ih[CI+1]][jh[j]][kh[k]];
            H.fnew[iH][jH][kH] += h.A.w[ih[1   ]][jh[j]][kh[k]] 
                                * h.phi[ih[0   ]][jh[j]][kh[k]];

            for(int i=1; i< CI; i++)
              D[i][j][k] -= (h.A.e[ih[i]][jh[j]][kh[k]]);

            for(int i=CI; i>1; i--)
              D[i][j][k] -= (h.A.w[ih[i]][jh[j]][kh[k]]);
          }

        /*------+
        |  n-s  |
        +------*/
        for(int i=1; i<=CI; i++)
          for(int k=1; k<=CK; k++) {
            H.fnew[iH][jH][kH] += h.A.n[ih[i]][jh[CJ  ]][kh[k]] 
                                * h.phi[ih[i]][jh[CJ+1]][kh[k]];
            H.fnew[iH][jH][kH] += h.A.s[ih[i]][jh[1   ]][kh[k]] 
                                * h.phi[ih[i]][jh[0   ]][kh[k]];
                                
            for(int j=1; j< CJ; j++)
              D[i][j][k] -= (h.A.n[ih[i]][jh[j]][kh[k]]);

            for(int j=CJ; j>1; j--)
              D[i][j][k] -= (h.A.s[ih[i]][jh[j]][kh[k]]);
          }

        /*------+
        |  t-b  |
        +------*/
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++) {
            H.fnew[iH][jH][kH] += h.A.t[ih[i]][jh[j]][kh[CK  ]]
                                * h.phi[ih[i]][jh[j]][kh[CK+1]];
            H.fnew[iH][jH][kH] += h.A.b[ih[i]][jh[j]][kh[1   ]]
                                * h.phi[ih[i]][jh[j]][kh[0   ]];
                
            for(int k=1; k< CK; k++)
              D[i][j][k] -= h.A.t[ih[i]][jh[j]][kh[k]];

            for(int k=CK; k> 1; k--)
              D[i][j][k] -= h.A.b[ih[i]][jh[j]][kh[k]];
          }

        /*---------------------------+ 
        |  finalize right hand side  |
        +---------------------------*/ 
        for(int i=1; i<=CI; i++)
          for(int j=1; j<=CJ; j++) 
            for(int k=1; k<=CK; k++)
              H.fnew[iH][jH][kH] -= D[i][j][k] * h.phi[ih[i]][jh[j]][kh[k]];
        }
      }
    }

#if 0 /* this does not appear to do anything useful ?! - Lubomir */
  /*----------------------------------------------------+ 
  |  r.h.s. normalization.                              |
  + - - - - - - - - - - - - - - - - - - - - - - - - - - +
  |  maybe this is quite normal. a.c. should solve      |
  |  coarse levels directly. if i do not do that, i     |
  |  get distorted r.h.s. on even coarser levels        |
  +----------------------------------------------------*/
  real sum = 0.0;
  for(int i=H.si(); i<=H.ei(); i++)
    for(int j=H.sj(); j<=H.ej(); j++)
      for(int k=H.sk(); k<=H.ek(); k++)
        sum += H.fnew[i][j][k];

  boil::cart.sum_real(& sum);

  int ncell = (H.ni()-2*boil::BW)*(H.nj()-2*boil::BW)*(H.nk()-2*boil::BW);
  boil::cart.sum_int( & ncell );

  sum /= (real)ncell;

  for(int i=H.si(); i<=H.ei(); i++)
    for(int j=H.sj(); j<=H.ej(); j++)
      for(int k=H.sk(); k<=H.ek(); k++)
        H.fnew[i][j][k] -= sum;
#endif

}
