#include "custom.h"

namespace boil {
  /******************************************************************************/
  void prolongate_vf_XZ(const Scalar & coarse, Scalar & fine, 
                        const VOF & concc, VOF & concf) {
  /***************************************************************************//**
   \brief Prolongate volume fraction geometrically from a coarse 2D grid 
          to a fine one.
      output: fine

     ---------     b = boil::BW-1
     | x | y |     coarse coords: i,j,k ; i^ = i-b, k^ = k-b
     ---------     fine coords:
     | q | z |         y - 2*i^   +b, j ,2*k^   +b
     ---------         x - 2*i^-1 +b, j ,2*k^   +b
                       z - 2*i^   +b, j ,2*k^-1 +b
                       q - 2*i^-1 +b, j ,2*k^-1 +b

  *******************************************************************************/

    for_vijk(coarse,i,j,k) {

      /* step 1: establish coordinate correspondence */
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

      /* step 2: prolongate */
      fine[i_x][j_x][k_x] = translate_v_coarse_to_fine(i_x,j_x,k_x,i,j,k,
                                                       coarse, fine, 
                                                       concc, concf);
      fine[i_y][j_y][k_y] = translate_v_coarse_to_fine(i_y,j_y,k_y,i,j,k,
                                                       coarse, fine, 
                                                       concc, concf);
      fine[i_z][j_z][k_z] = translate_v_coarse_to_fine(i_z,j_z,k_z,i,j,k,
                                                       coarse, fine, 
                                                       concc, concf);
      fine[i_q][j_q][k_q] = translate_v_coarse_to_fine(i_q,j_q,k_q,i,j,k,
                                                       coarse, fine, 
                                                       concc, concf);
    }

    fine.bnd_update();
    fine.exchange_all();

    concf.nx.bnd_update();
    concf.ny.bnd_update();
    concf.nz.bnd_update();
    concf.nalpha.bnd_update();
    concf.nx.exchange_all();
    concf.ny.exchange_all();
    concf.nz.exchange_all();
    concf.nalpha.exchange_all();

    return;
  }

  inline real translate_v_coarse_to_fine(const int i_f, const int j_f, const int k_f,
                                         const int i_c, const int j_c, const int k_c,
                                         const Scalar & coarse, const Scalar & fine,
                                         const VOF & concc, VOF & concf) {
    return concc.translate_v(i_c,j_c,k_c,
                             fine.xc(i_f)-coarse.xc(i_c),
                             fine.yc(j_f)-coarse.yc(j_c),
                             fine.zc(k_f)-coarse.zc(k_c),
                             fine.dxc(i_f)/coarse.dxc(i_c),
                             fine.dyc(j_f)/coarse.dyc(j_c),
                             fine.dzc(k_f)/coarse.dzc(k_c),
                             concf.nx[i_f][j_f][k_f],
                             concf.ny[i_f][j_f][k_f],
                             concf.nz[i_f][j_f][k_f],
                             concf.nalpha[i_f][j_f][k_f],
                             coarse);
  }

}
