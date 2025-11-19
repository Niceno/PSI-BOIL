#include "custom.h"
  
namespace boil {
  /******************************************************************************/
  void output_wall_heat_transfer_xz(const CommonHeatTransfer & cht,
                                    std::ostream & otp, const int NX) {
  /***************************************************************************//**
    \brief Output wall temperature and heat flux in z-direction, assuming solid
           wall at z = 0 to an output stream. 
           NX = total number of cells in x. 
           Transformation i->glob_i-boil::BW+1 is used,
           result: closed interval [1,NX]
  *******************************************************************************/
    for(int I(1); I<=NX; ++I) {
      real tprsum0(-boil::unreal), tprsum1(-boil::unreal),
           bndtprsum0(-boil::unreal), bndtprsum1(-boil::unreal),
           hf(-boil::unreal), xpos(-boil::unreal);
      for_vi(cht.tmp(),i) {
        if(cht.tmp().domain()->global_I(i)==(I+boil::BW-1)) {
          for_vjk(cht.tmp(),j,k) {
            if(cht.tmp().zc(k)<0&&cht.tmp().zc(k+1)>0) {
              tprsum0 = cht.tmp()[i][j][k];
              bndtprsum0 = cht.node_tmp_sol()[Comp::w()][i][j][k+1];
              bndtprsum1 = cht.node_tmp_flu()[Comp::w()][i][j][k+1];
              tprsum1 = cht.tmp()[i][j][k+1];
              xpos = cht.tmp().xc(i);
              real lambdas = cht.lambda(i,j,k);
              real resist_s = 0.5*cht.tmp().dzc(k)/lambdas;
              hf = (tprsum0-bndtprsum0)/resist_s + cht.dirac_wall_source(i,j,k+1);
            }
          } /* vjk */
        } /* I is local i */
      } /* local i */
      boil::cart.max_real(&tprsum0);
      boil::cart.max_real(&bndtprsum0);
      boil::cart.max_real(&tprsum1);
      boil::cart.max_real(&bndtprsum1);
      boil::cart.max_real(&hf);
      boil::cart.max_real(&xpos);
      if(!boil::cart.iam()) {
        otp<<I<<" "<<xpos<<" "
              <<tprsum0<<" "<<bndtprsum0<<" "
              <<bndtprsum1<<" "<<tprsum1<<" "
              <<hf<<boil::endl;
      }
    } /* global I */

    return;
  }

  /******************************************************************************/
  void output_wall_heat_transfer_xz(const Scalar & tpr,
                                    const Topology & topo,
                                    const PhaseChange4 & pc,
                                    std::ostream & otp, const int NX) {
  /***************************************************************************//**
    \brief Output wall temperature and heat flux in z-direction, assuming solid
           wall at z = 0 to an output stream, no conjugate HT. 
           NX = total number of cells in x. 
           Transformation i->glob_i-boil::BW+1 is used,
           result: closed interval [1,NX]
  *******************************************************************************/

    for(int I(1); I<=NX; ++I) {
      real tprsum(-boil::unreal), bndtprsum(-boil::unreal),
           hf1(-boil::unreal), hf2(-boil::unreal), xpos(-boil::unreal);
      for_vi(tpr,i) {
        if(tpr.domain()->global_I(i)==(I+boil::BW-1)) {
          for_vjk(tpr,j,k) {
            if(tpr.zc(k-1)<boil::atto&&tpr.zc(k)>0) {
              tprsum = tpr[i][j][k];
              bndtprsum = tpr[i][j][k-1];
              std::vector<real> tvs, tls;
              std::vector<real> tvf, tlf;
              pc.request_flux(i,j,k-1,tvs,tls);
              pc.request_flux(i,j,k,tvf,tlf);
              xpos = tpr.xc(i);
              if(topo.above_interface(i,j,k)) {  
                if(topo.interface(Sign::neg(),Comp::k(),i,j,k)) {
                  hf1 = -tvs[2];
                  hf2 = -tvf[2];
                } else {
                  hf1 = -tls[2];
                  hf2 = -tlf[2];
                }
              } else {
                if(topo.interface(Sign::neg(),Comp::k(),i,j,k)) {
                  hf1 = -tls[2];
                  hf2 = -tlf[2];
                } else {
                  hf1 = -tvs[2];
                  hf2 = -tvf[2];
                }
              } /* above interface */
            } /* next to wall */
          } /* vjk */
        } /* I is local i */
      } /* local i */
      boil::cart.max_real(&tprsum);
      boil::cart.max_real(&bndtprsum);
      boil::cart.max_real(&hf1);
      boil::cart.max_real(&hf2);
      boil::cart.max_real(&xpos);
      if(!boil::cart.iam()) {
        otp<<I<<" "<<xpos<<" "<<bndtprsum<<" "<<tprsum<<" "<<hf1<<" "<<hf2<<boil::endl;
      }
    } /* global I */

    return;
  }
}
