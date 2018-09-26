#include "preconditioner.h"

/******************************************************************************/
void Diagonal :: form(const Matrix & A, const Scalar & x) {
/*--------------------------------+
|  form preconditioning matrix M  |
+--------------------------------*/

  Mc = x.shape(); Mc=0.0;

  for_vijk(Mc,i,j,k)
    Mc[i][j][k] = 1.0/A.c[i][j][k];
}  

/******************************************************************************/
void IncompleteCholesky0 :: form(const Matrix & A, const Scalar & x) {
/*--------------------------------+
|  form preconditioning matrix M  |
+--------------------------------*/

  Mc = x.shape(); Mc=0.0;
  Mw = x.shape(); Mw=0.0;
  Ms = x.shape(); Ms=0.0;
  Mb = x.shape(); Mb=0.0;

  for_vijk(Mc,i,j,k) {
    real sum = A.c[i][j][k];
    sum = sum - Mw[i][j][k] * Mw[i][j][k];
    sum = sum - Ms[i][j][k] * Ms[i][j][k];
    sum = sum - Mb[i][j][k] * Mb[i][j][k];
    sum = sqrt(sum);
    Mc[i][j][k] = sum;

    Mw[i+1][j][k] = - A.w[i+1][j][k]/Mc[i][j][k];
    Ms[i][j+1][k] = - A.s[i][j+1][k]/Mc[i][j][k];
    Mb[i][j][k+1] = - A.b[i][j][k+1]/Mc[i][j][k];
  }

  // one must not exchange like:
  //Mw.exchange();
  //Ms.exchange();
  //Mb.exchange();

}  

/******************************************************************************/
void IncompleteCholesky2 :: form(const Matrix & A, const Scalar & x) {
/*--------------------------------+
|  form preconditioning matrix M  |
+--------------------------------*/
 
  Mc  = x.shape(); Mc =0.0;
  Mw  = x.shape(); Mw =0.0;
  Ms  = x.shape(); Ms =0.0;
  Mb  = x.shape(); Mb =0.0;
  Mtw = x.shape(); Mtw=0.0;
  Mts = x.shape(); Mts=0.0;

  for_vijk(Mc,i,j,k) {
    real sum = A.c[i][j][k];
    sum = sum - Mw [i][j][k] * Mw [i][j][k];
    sum = sum - Mtw[i][j][k] * Mtw[i][j][k];
    sum = sum - Ms [i][j][k] * Ms [i][j][k];
    sum = sum - Mts[i][j][k] * Mts[i][j][k];
    sum = sum - Mb [i][j][k] * Mb [i][j][k];
    sum = sqrt(sum);
    Mc[i][j][k] = sum;

    Mw [i+1][j]  [k]   = - A.w[i+1][j][k]/Mc[i][j][k];

    Mtw[i+1][j]  [k-1] = - Mw[i+1][j][k-1]*Mb[i][j][k]/Mc[i][j][k];

    Ms [i]  [j+1][k]   = - A.s[i][j+1][k]/Mc[i][j][k];

    Mts[i]  [j+1][k-1] = - Ms[i][j+1][k-1]*Mb[i][j][k]/Mc[i][j][k];

    Mb [i]  [j]  [k+1] = (- A.b[i][j][k+1]
                          - Mw[i][j][k+1]*Mtw[i][j][k]
                          - Ms[i][j][k+1]*Mts[i][j][k])/Mc[i][j][k];
  }

  // one must not exchange like:
  //Mw.exchange();
  //Mtw.exchange();
  //Ms.exchange();
  //Mts.exchange();
  //Mb.exchange();

}  

/******************************************************************************/
void IncompleteCholesky3 :: form(const Matrix & A, const Scalar & x) {
/*--------------------------------+
|  form preconditioning matrix M  |
+--------------------------------*/
 
  Mc  = x.shape(); Mc =0.0;
  Mw  = x.shape(); Mw =0.0;
  Ms  = x.shape(); Ms =0.0;
  Mb  = x.shape(); Mb =0.0;
  Mtw = x.shape(); Mtw=0.0;
  Mts = x.shape(); Mts=0.0;
  Mnw = x.shape(); Mnw=0.0;

  for_vijk(Mc,i,j,k) {
    real sum = A.c[i][j][k];
    sum = sum - Mw [i][j][k] * Mw [i][j][k];
    sum = sum - Mtw[i][j][k] * Mtw[i][j][k];
    sum = sum - Mnw[i][j][k] * Mnw[i][j][k];
    sum = sum - Ms [i][j][k] * Ms [i][j][k];
    sum = sum - Mts[i][j][k] * Mts[i][j][k];
    sum = sum - Mb [i][j][k] * Mb [i][j][k];
    sum = sqrt(sum);
    Mc[i][j][k] = sum;

    Mw [i+1][j]  [k]   = - A.w[i+1][j][k]/Mc[i][j][k];

    Mtw[i+1][j]  [k-1] = - Mw[i+1][j][k-1]*Mb[i][j][k]/Mc[i][j][k];

    Mnw[i+1][j-1][k]   = (- Mw [i+1][j-1][k]*Ms [i][j][k] 
                          - Mtw[i+1][j-1][k]*Mts[i][j][k])/Mc[i][j][k];

    Ms [i]  [j+1][k]   = (- A.s[i][j+1][k]
                          - Mw [i][j+1][k]*Mnw[i][j][k]
                          - Mnw[i][j+1][k]*Ms [i][j][k])/Mc[i][j][k];

    Mts[i]  [j+1][k-1] = (- Ms [i][j+1][k-1]*Mb [i][j][k]
                          - Mtw[i][j+1][k-1]*Mnw[i][j][k])/Mc[i][j][k];

    Mb [i]  [j]  [k+1] = (- A.b[i][j][k+1]
                          - Mw[i][j][k+1]*Mtw[i][j][k]
                          - Ms[i][j][k+1]*Mts[i][j][k])/Mc[i][j][k];
  }

  // one must not exchange like:
  //Mw.exchange();
  //Mtw.exchange();
  //Mnw.exchange();
  //Ms.exchange();
  //Ms.exchange();
  //Mts.exchange();
  //Mb.exchange();

}  

/*-----------------------------------------------------------------------------+
 '$Id: preconditioner_form.cpp,v 1.16 2016/02/12 13:49:55 sato Exp $'/
+-----------------------------------------------------------------------------*/
