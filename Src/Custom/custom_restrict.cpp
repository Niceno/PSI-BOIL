#include "custom.h"

namespace boil {
  /******************************************************************************/
  void restrictXZ(const Scalar & fine, Scalar & coarse) {
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
  void restrictXZ_area(const Scalar & fine, Scalar & coarse) {
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

      if(fine.domain()->is_cartesian()) {
        coefx = fine.dSy(i_x,j_x,k_x);
        coefy = fine.dSy(i_y,j_y,k_y);
        coefz = fine.dSy(i_z,j_z,k_z);
        coefq = fine.dSy(i_q,j_q,k_q);

        coefcoarse = coarse.dSy(i,j,k);
      } else {
        coefx = fine.domain()->dSy_cartesian(i_x,j_x,k_x);
        coefy = fine.domain()->dSy_cartesian(i_y,j_y,k_y);
        coefz = fine.domain()->dSy_cartesian(i_z,j_z,k_z);
        coefq = fine.domain()->dSy_cartesian(i_q,j_q,k_q);

        coefcoarse = coarse.domain()->dSy_cartesian(i,j,k);
      }

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
  void restrictXZ_vector_simple(const Vector & fine, Vector & coarse,
                                const Scalar & fs, const Scalar & cs) {
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

      real coefx = fs.dSx(Sign::neg(),i_x,j_x,k_x);
      real coefq = fs.dSx(Sign::neg(),i_q,j_q,k_q);
      real coefy = fs.dSx(Sign::pos(),i_y,j_y,k_y);
      real coefz = fs.dSx(Sign::pos(),i_z,j_z,k_z);

      real coefc = cs.dSx(Sign::neg(),i,j,k);

      if(coefc<boil::atto) {
        coefx = coefq = coefc = 1.0;
      }

      coarse[m][i  ][j][k] = fine[m][i_x][j_x][k_x]*coefx
                           + fine[m][i_q][j_q][k_q]*coefq;

      coarse[m][i  ][j][k] /= coefc;

      coefc = cs.dSx(Sign::pos(),i,j,k);

      if(coefc<boil::atto) {
        coefy = coefz = coefc = 1.0;
      }

      coarse[m][i+1][j][k] = fine[m][i_y][j_y][k_y]*coefy
                           + fine[m][i_z][j_z][k_z]*coefz;

      coarse[m][i+1][j][k] /= coefc;

      m = Comp::w();

      coefx = fs.dSz(Sign::pos(),i_x,j_x,k_x);
      coefq = fs.dSz(Sign::neg(),i_q,j_q,k_q);
      coefy = fs.dSz(Sign::pos(),i_y,j_y,k_y);
      coefz = fs.dSz(Sign::neg(),i_z,j_z,k_z);

      coefc = cs.dSz(Sign::neg(),i,j,k);

      if(coefc<boil::atto) {
        coefq = coefz = coefc = 1.0;
      }

      coarse[m][i][j][k  ] = fine[m][i_q][j_q][k_q]*coefq
                           + fine[m][i_z][j_z][k_z]*coefz;
 
      coarse[m][i][j][k  ] /= coefc;

      coefc = cs.dSz(Sign::pos(),i,j,k);

      if(coefc<boil::atto) {
        coefx = coefy = coefc = 1.0;
      }

      coarse[m][i][j][k+1] = fine[m][i_x][j_x][k_x]*coefx
                           + fine[m][i_y][j_y][k_y]*coefy;

      coarse[m][i][j][k+1] /= coefc;

    }
    coarse.exchange_all();

    return;
  }


}
