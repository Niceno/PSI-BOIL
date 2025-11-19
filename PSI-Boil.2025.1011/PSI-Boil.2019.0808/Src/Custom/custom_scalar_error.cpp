#include "custom.h"

namespace boil {
  /******************************************************************************/
  real l2_scalar_error(const Scalar & sca, const Scalar & scb) {
  /***************************************************************************//**
   \brief Calculate the L2 difference between two scalar fields.
      output: L2
  *******************************************************************************/

    real l2(-0.0);
    int cnt(0);
    for_vijk(sca,i,j,k) {
      //if(dom->ibody().off(i,j,k)) continue;
      //real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
      //                 -std::max(0.0,std::min(1.0,scb[i][j][k])));
      real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
      l2 += sdiff*sdiff;
      cnt++;
    }
    boil::cart.sum_real(&l2);
    boil::cart.sum_int(&cnt);

    return l2;
  }

  /******************************************************************************/
  real l1_scalar_error(const Scalar & sca, const Scalar & scb) {
  /***************************************************************************//**
   \brief Calculate the L1 difference between two scalar fields.
      output: L1
  *******************************************************************************/

    real l1(-0.0);
    int cnt(0);
    for_vijk(sca,i,j,k) {
      real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
      l1 += sdiff;
      cnt++;
    }
    boil::cart.sum_real(&l1);
    boil::cart.sum_int(&cnt);

    return l1;
  }

  /******************************************************************************/
  real li_scalar_error(const Scalar & sca, const Scalar & scb) {
  /***************************************************************************//**
   \brief Calculate the Li difference between two scalar fields.
      output: Li
  *******************************************************************************/

    real li(-0.0);
    int cnt(0);
    for_vijk(sca,i,j,k) {
      //if(dom->ibody().off(i,j,k)) continue;
      //real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
      //                 -std::max(0.0,std::min(1.0,scb[i][j][k])));
      real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
      if(sdiff>li)
        li = sdiff;
    }
    boil::cart.max_real(&li);

    return li;
  }

  /******************************************************************************/
  real l2_scalar_error_vol(const Scalar & sca, const Scalar & scb) {
  /***************************************************************************//**
   \brief Calculate the L2 difference between two scalar fields.
      output: L2
  *******************************************************************************/

    real l2(-0.0);
    real vol(0.);
    for_vijk(sca,i,j,k) {
      //if(dom->ibody().off(i,j,k)) continue;
      //real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
      //                 -std::max(0.0,std::min(1.0,scb[i][j][k])));
      real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
      l2 += sdiff*sdiff*sca.dV(i,j,k);
      vol += sca.dV(i,j,k);
    }
    boil::cart.sum_real(&l2);
    boil::cart.sum_real(&vol);

    return sqrt(l2/vol);
  }

  /******************************************************************************/
  real l1_scalar_error_vol(const Scalar & sca, const Scalar & scb) {
  /***************************************************************************//**
   \brief Calculate the L1 difference between two scalar fields.
      output: L1
  *******************************************************************************/

    real l1(-0.0);
    real vol(0.);
    for_vijk(sca,i,j,k) {
      real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
      l1 += sdiff*sca.dV(i,j,k);
      vol += sca.dV(i,j,k);
    }
    boil::cart.sum_real(&l1);
    boil::cart.sum_real(&vol);

    return l1/vol;
  }
}
