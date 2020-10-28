#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  void output_profile_xz(const Scalar & c, std::ostream & otp,
                         const Range<int> RZ, const Range<int> RX,
                         const real subtract) {
  /***************************************************************************//**
    \brief Integrate 2D scalar in the x-z plane for all layers in z-direction and
           output the result to an output stream. 
           RZ = range of cells to be considered in Z. 
           RX = range of cells to be considered in X (default is all).
           Transformation k->glob_k-boil::BW+1 is used,
           result: closed interval based on the range
  *******************************************************************************/

    bool ranged(false);
    if(RX.exists()) {
      ranged = true;
    }

    for(int K(RZ.first()); K<=RZ.last(); ++K) {
      real xsum(0.0),zpos(0.0);
      for_vk(c,k) {
        if(c.domain()->global_K(k)==(K+boil::BW-1)) {
          zpos = c.zc(k);
          for_vi(c,i) {
            if(!ranged||RX.contains(c.domain()->global_I(i))) {
              for_vj(c,j) { 
                if(c.domain()->ibody().on(i,j,k)) {
                  xsum += c[i][j][k]*c.dxc(i);
                }
              } /* vj = should be just one layer */
            } /* i is in range */
          } /* vi */
        } /* K is local k */
      } /* vk*/
      boil::cart.max_real(&zpos);
      boil::cart.sum_real(&xsum);
      if(subtract>0.)
        xsum = subtract - xsum;
      if(!boil::cart.iam()) {
        otp<<K<<" "<<zpos<<" "<<xsum<<boil::endl;
      }
    } /* all K */

    return;
  }

  /******************************************************************************/
  void output_profile_zx(const Scalar & c, std::ostream & otp,
                         const Range<int> RX, const Range<int> RZ,
                         const real subtract) {
  /***************************************************************************//**
    \brief Integrate 2D scalar in the x-z plane for all layers in x-direction and
           output the result to an output stream. 
           RX = range of cells to be considered in X. 
           RZ = range of cells to be considered in Z (default is all).
           Transformation i->glob_i-boil::BW+1 is used,
           result: closed interval based on the range
  *******************************************************************************/

    bool ranged(false);
    if(RX.exists()) {
      ranged = true;
    }

    for(int I(RX.first()); I<=RX.last(); ++I) {
      real zsum(0.0),xpos(0.0);
      for_vi(c,i) {
        if(c.domain()->global_I(i)==(I+boil::BW-1)) {
          xpos = c.xc(i);
          for_vk(c,k) {
            if(!ranged||RZ.contains(c.domain()->global_K(k))) {
              for_vj(c,j) {
                if(c.domain()->ibody().on(i,j,k)) {
                  zsum += c[i][j][k]*c.dzc(k);
                }
              } /* vj = should be just one layer */
            } /* k is in range */
          } /* vk */
        } /* I is local i */
      } /* vi*/
      boil::cart.max_real(&xpos);
      boil::cart.sum_real(&zsum);
      if(subtract>0.)
        zsum = subtract - zsum;
      if(!boil::cart.iam()) {
        otp<<I<<" "<<xpos<<" "<<zsum<<boil::endl;
      }
    } /* all I */

    return;
  }
}
