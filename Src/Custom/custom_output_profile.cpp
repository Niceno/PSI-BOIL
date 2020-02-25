#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  void output_profile_xz(const Scalar & c, std::ostream & otp, const int NZ) {
  /***************************************************************************//**
    \brief Integrate 2D scalar in the x-z plane for all layers in z-direction and
           output the result to an output stream. 
           NZ = total number of cells in z. 
           Transformation k->glob_k-boil::BW+1 is used,
           result: closed interval [1,NZ]
  *******************************************************************************/
    for(int K(1); K<=NZ; ++K) {
      real xsum(0.0),zpos(0.0);
      for_vk(c,k) {
        if(c.domain()->global_K(k)==(K+boil::BW-1)) {
          zpos = c.zc(k);
          for_vij(c,i,j) { 
            if(c.domain()->ibody().on(i,j,k)) {
              xsum += c[i][j][k]*c.dxc(i);
            }
          } /* vij */
        } /* K is local k */
      } /* vk*/
      boil::cart.max_real(&zpos);
      boil::cart.sum_real(&xsum);
      if(!boil::cart.iam()) {
        otp<<K<<" "<<zpos<<" "<<xsum<<boil::endl;
      }
    } /* all K */

    return;
  }
}
