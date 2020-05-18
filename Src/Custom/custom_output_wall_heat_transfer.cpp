#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  void output_wall_heat_transfer_xz(const Scalar & tpr,
                                    const Vector & bndtpr,
                                    const real lambdas,
                                    std::ostream & otp, const int NX) {
  /***************************************************************************//**
    \brief Output wall temperature and heat flux in z-direction, assuming solid
           wall at z = 0 to an output stream. 
           NX = total number of cells in x. 
           Transformation i->glob_i-boil::BW+1 is used,
           result: closed interval [1,NX]
  *******************************************************************************/
    for(int I(1); I<=NX; ++I) {
      real tprsum(-boil::unreal), bndtprsum(-boil::unreal),
           hf(-boil::unreal), xpos(-boil::unreal);
      for_vi(tpr,i) {
        if(tpr.domain()->global_I(i)==(I+boil::BW-1)) {
          for_vjk(tpr,j,k) {
            if(tpr.zc(k)<0&&tpr.zc(k+1)>0) {
              tprsum = tpr[i][j][k];
              bndtprsum = bndtpr[Comp::w()][i][j][k+1];
              xpos = tpr.xc(i);
              hf = lambdas*(tprsum-bndtprsum)/(0.5*tpr.dzc(k));
            }
          } /* vjk */
        } /* I is local i */
      } /* local i */
      boil::cart.max_real(&tprsum);
      boil::cart.max_real(&bndtprsum);
      boil::cart.max_real(&hf);
      boil::cart.max_real(&xpos);
      if(!boil::cart.iam()) {
        otp<<I<<" "<<xpos<<" "<<tprsum<<" "<<bndtprsum<<" "<<hf<<boil::endl;
      }
    } /* global I */

    return;
  }
}