#include "preconditioner.h"

/******************************************************************************/
void Diagonal :: solve(Scalar & z, const Scalar & r) {
/*-------------------------------+
|  solve preconditioning sjstem  |
+-------------------------------*/

  for_vijk(Mc,i,j,k)
    z[i][j][k] = r[i][j][k] * Mc[i][j][k];
}  

/******************************************************************************/
void IncompleteCholesky0 :: solve(Scalar & z, const Scalar & r) {
/*-------------------------------+
|  solve preconditioning sjstem  |
+-------------------------------*/

  for(int i=Mc.si(); i<=Mc.ei(); i++)
    for(int j=Mc.sj(); j<=Mc.ej(); j++)
      for(int k=Mc.sk(); k<=Mc.ek(); k++) {
        real sum = r[i][j][k];
        sum = sum - Mw[i][j][k]*z[i-1][j][k];
        sum = sum - Ms[i][j][k]*z[i][j-1][k];
        sum = sum - Mb[i][j][k]*z[i][j][k-1];
        z[i][j][k] = sum/Mc[i][j][k];
      }

  for(int i=Mc.ei(); i>=Mc.si(); i--)
    for(int j=Mc.ej(); j>=Mc.sj(); j--)
      for(int k=Mc.ek(); k>=Mc.sk(); k--) {
        real sum = z[i][j][k];
        sum = sum - Mw[i+1][j][k]*z[i+1][j][k];
        sum = sum - Ms[i][j+1][k]*z[i][j+1][k];
        sum = sum - Mb[i][j][k+1]*z[i][j][k+1];
        z[i][j][k] = sum/Mc[i][j][k];
      }
}  

/******************************************************************************/
void IncompleteCholesky2 :: solve(Scalar & z, const Scalar & r) {
/*-------------------------------+
|  solve preconditioning sjstem  |
+-------------------------------*/

  for(int i=Mc.si(); i<=Mc.ei(); i++)
    for(int j=Mc.sj(); j<=Mc.ej(); j++)
      for(int k=Mc.sk(); k<=Mc.ek(); k++) {
        real sum = r[i][j][k];

        sum = sum - Mw [i][j][k] * z[i-1][j]  [k];
        sum = sum - Mtw[i][j][k] * z[i-1][j]  [k+1];
        sum = sum - Ms [i][j][k] * z[i]  [j-1][k];
        sum = sum - Mts[i][j][k] * z[i]  [j-1][k+1];
        sum = sum - Mb [i][j][k] * z[i]  [j]  [k-1];

        z[i][j][k] = sum/Mc[i][j][k];
      }

  for(int i=Mc.ei(); i>=Mc.si(); i--)
    for(int j=Mc.ej(); j>=Mc.sj(); j--)
      for(int k=Mc.ek(); k>=Mc.sk(); k--) {
        real sum = z[i][j][k];

        sum = sum - Mw [i+1][j]  [k]   * z[i+1][j]  [k];
        sum = sum - Mtw[i+1][j]  [k-1] * z[i+1][j]  [k-1];
        sum = sum - Ms [i]  [j+1][k]   * z[i]  [j+1][k];
        sum = sum - Mts[i]  [j+1][k-1] * z[i]  [j+1][k-1];
        sum = sum - Mb [i]  [j]  [k+1] * z[i]  [j]  [k+1];

        z[i][j][k] = sum/Mc[i][j][k];
      }
}  

/******************************************************************************/
void IncompleteCholesky3 :: solve(Scalar & z, const Scalar & r) {
/*-------------------------------+
|  solve preconditioning sjstem  |
+-------------------------------*/

  for(int i=Mc.si(); i<=Mc.ei(); i++)
    for(int j=Mc.sj(); j<=Mc.ej(); j++)
      for(int k=Mc.sk(); k<=Mc.ek(); k++) {
        real sum = r[i][j][k];

        sum = sum - Mw [i][j][k] * z[i-1][j]  [k];
        sum = sum - Mtw[i][j][k] * z[i-1][j]  [k+1];
        sum = sum - Mnw[i][j][k] * z[i-1][j+1][k];  
        sum = sum - Ms [i][j][k] * z[i]  [j-1][k];
        sum = sum - Mts[i][j][k] * z[i]  [j-1][k+1];
        sum = sum - Mb [i][j][k] * z[i]  [j]  [k-1];

        z[i][j][k] = sum/Mc[i][j][k];
      }

  for(int i=Mc.ei(); i>=Mc.si(); i--)
    for(int j=Mc.ej(); j>=Mc.sj(); j--)
      for(int k=Mc.ek(); k>=Mc.sk(); k--) {
        real sum = z[i][j][k];

        sum = sum - Mw [i+1][j]  [k]   * z[i+1][j]  [k];
        sum = sum - Mtw[i+1][j]  [k-1] * z[i+1][j]  [k-1];
        sum = sum - Mnw[i+1][j-1][k]   * z[i+1][j-1][k];  
        sum = sum - Ms [i]  [j+1][k]   * z[i]  [j+1][k];
        sum = sum - Mts[i]  [j+1][k-1] * z[i]  [j+1][k-1];
        sum = sum - Mb [i]  [j]  [k+1] * z[i]  [j]  [k+1];

        z[i][j][k] = sum/Mc[i][j][k];
      }
}  

/*-----------------------------------------------------------------------------+
 '$Id: preconditioner_solve.cpp,v 1.10 2011/09/29 03:12:57 niceno Exp $'/
+-----------------------------------------------------------------------------*/
