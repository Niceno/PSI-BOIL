#include "twolevelvector.h"

/* only 2D X-Z version implemented */
/******************************************************************************/
void TwoLevelVector::restrict_area_XZ(const Scalar & cs) {
/******************************************************************************/

  for_vijk(cs,i,j,k) {
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
    
    Comp m = Comp::u();

    real coefx = fine.dSx(m,i_x,j_x,k_x);
    real coefq = fine.dSx(m,i_q,j_q,k_q);
    real coefy = fine.dSx(m,i_y+1,j_y,k_y);
    real coefz = fine.dSx(m,i_z+1,j_z,k_z);

    real coefc = coarse.dSx(m,i,j,k);

    if(coefc<boil::atto) {
      coefx = coefq = coefc = 1.0;
    }

    coarse[m][i  ][j][k] = fine[m][i_x][j_x][k_x]*coefx
                         + fine[m][i_q][j_q][k_q]*coefq;

    coarse[m][i  ][j][k] /= coefc;

    coefc = coarse.dSx(m,i+1,j,k);

    if(coefc<boil::atto) {
      coefy = coefz = coefc = 1.0;
    }

    coarse[m][i+1][j][k] = fine[m][i_y+1][j_y][k_y]*coefy
                         + fine[m][i_z+1][j_z][k_z]*coefz;

    coarse[m][i+1][j][k] /= coefc;

    m = Comp::w();

    coefx = fine.dSz(m,i_x,j_x,k_x);
    coefq = fine.dSz(m,i_q,j_q,k_q);
    coefy = fine.dSz(m,i_y,j_y,k_y+1);
    coefz = fine.dSz(m,i_z,j_z,k_z+1);

    coefc = coarse.dSz(m,i,j,k);

    if(coefc<boil::atto) {
      coefq = coefz = coefc = 1.0;
    }

    coarse[m][i][j][k  ] = fine[m][i_q][j_q][k_q]*coefq
                         + fine[m][i_z][j_z][k_z]*coefz;

    coarse[m][i][j][k  ] /= coefc;

    coefc = coarse.dSz(m,i,j,k+1);

    if(coefc<boil::atto) {
      coefx = coefy = coefc = 1.0;
    }

    coarse[m][i][j][k+1] = fine[m][i_x][j_x][k_x+1]*coefx
                         + fine[m][i_y][j_y][k_y+1]*coefy;

    coarse[m][i][j][k+1] /= coefc;

  }
  coarse.exchange_all();

  return;
}
