#include "twolevelscalar.h"

/* only 2D X-Z versions implemented */

/******************************************************************************/
void TwoLevelScalar::restrict_volume_XZ() {
/***************************************************************************//**
 \brief Restrict scalar on a 2D fine grid to a coarse one.
    output: coarse

   ---------     b = boil::BW-1
   | x | y |     coarse coords: i,j,k ; i^ = i-b, k^ = k-b
   ---------     fine coords:
   | q | z |         y - 2*i^   +b, j ,2*k^   +b
   ---------         x - 2*i^-1 +b, j ,2*k^   +b
                     z - 2*i^   +b, j ,2*k^-1 +b
                     q - 2*i^-1 +b, j ,2*k^-1 +b

*******************************************************************************/

  for_vijk(coarse,i,j,k) {
    /* establish coordinate correspondence */
    const int b = boil::BW-1;
    const int ihat = i-b;
    const int khat = k-b;

    const int i_y = 2*ihat + b;
    const int j_y = j;
    const int k_y = 2*khat + b;

    const int i_x = 2*ihat-1 + b;
    const int j_x = j;
    const int k_x = 2*khat + b;

    const int i_z = 2*ihat + b;
    const int j_z = j;
    const int k_z = 2*khat-1 + b;

    const int i_q = 2*ihat-1 + b;
    const int j_q = j;
    const int k_q = 2*khat-1 + b;

    coarse[i][j][k] = fine[i_x][j_x][k_x]*fine.dV(i_x,j_x,k_x)
                    + fine[i_y][j_y][k_y]*fine.dV(i_y,j_y,k_y)
                    + fine[i_z][j_z][k_z]*fine.dV(i_z,j_z,k_z)
                    + fine[i_q][j_q][k_q]*fine.dV(i_q,j_q,k_q);

    coarse[i][j][k] /= coarse.dV(i,j,k);
  }
  coarse.exchange_all();

  return;
}

/******************************************************************************/
void TwoLevelScalar::restrict_area_XZ() {
/******************************************************************************/

  for_vijk(coarse,i,j,k) {
    /* establish coordinate correspondence */
    const int b = boil::BW-1;
    const int ihat = i-b;
    const int khat = k-b;

    const int i_y = 2*ihat + b;
    const int j_y = j;
    const int k_y = 2*khat + b;

    const int i_x = 2*ihat-1 + b;
    const int j_x = j;
    const int k_x = 2*khat + b;

    const int i_z = 2*ihat + b;
    const int j_z = j;
    const int k_z = 2*khat-1 + b;

    const int i_q = 2*ihat-1 + b;
    const int j_q = j;
    const int k_q = 2*khat-1 + b;

    real coefx, coefy, coefz, coefq;
    real coefcoarse;

    /* protection against axisymmetry */
#if 0
    coefx = fine.dSy(i_x,j_x,k_x);
    coefy = fine.dSy(i_y,j_y,k_y);
    coefz = fine.dSy(i_z,j_z,k_z);
    coefq = fine.dSy(i_q,j_q,k_q);

    coefcoarse = coarse.dSy(i,j,k);
#else
    coefx = fine.domain()->dSy_cartesian(i_x,j_x,k_x);
    coefy = fine.domain()->dSy_cartesian(i_y,j_y,k_y);
    coefz = fine.domain()->dSy_cartesian(i_z,j_z,k_z);
    coefq = fine.domain()->dSy_cartesian(i_q,j_q,k_q);

    coefcoarse = coarse.domain()->dSy_cartesian(i,j,k);
#endif

    coarse[i][j][k] = fine[i_x][j_x][k_x]*coefx
                    + fine[i_y][j_y][k_y]*coefy
                    + fine[i_z][j_z][k_z]*coefz
                    + fine[i_q][j_q][k_q]*coefq;

    coarse[i][j][k] /= coefcoarse;

  }
  coarse.exchange_all();

  return;
}

/******************************************************************************/
void TwoLevelScalar::restrict_sum_XZ() {
/******************************************************************************/

  for_vijk(coarse,i,j,k) {
    /* establish coordinate correspondence */
    const int b = boil::BW-1;
    const int ihat = i-b;
    const int khat = k-b;

    const int i_y = 2*ihat + b;
    const int j_y = j;
    const int k_y = 2*khat + b;

    const int i_x = 2*ihat-1 + b;
    const int j_x = j;
    const int k_x = 2*khat + b;

    const int i_z = 2*ihat + b;
    const int j_z = j;
    const int k_z = 2*khat-1 + b;

    const int i_q = 2*ihat-1 + b;
    const int j_q = j;
    const int k_q = 2*khat-1 + b;

    coarse[i][j][k] = fine[i_x][j_x][k_x]
                    + fine[i_y][j_y][k_y]
                    + fine[i_z][j_z][k_z]
                    + fine[i_q][j_q][k_q];
  }
  coarse.exchange_all();

  return;
}
