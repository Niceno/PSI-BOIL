#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::m(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate M, usually in unit kg/m2s.
*         M = (qflux_liquid + qflux_vapor) / latent
*******************************************************************************/

#if 0 /* obtaining phase change rate from the underlying solid */
  for_ijk(i,j,k){
    if(Interface(i,j,k)){
      real lv = lambdav;
      real ll = lambdal;
      if (diff_eddy) {
        lv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        ll += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      /* z direction */
      Comp m = Comp::w();

      /* bottom is in wall & there is interface in bottom */
      if(dom->ibody().off(i,j,k-1)) {
        real tmp_w = bndtpr[m][i][j][k];
        real clrc = clr[i][j][k];
        real len_s = phi.dzb(k) - 0.5*phi.dzc(k);
        real lam_s = solid()->lambda(i,j,k-1);
        real tmp_s = tpr[i][j][k-1];

        if(Interface(-1,m,i,j,k)) {
          real lam_f;
          if(clrc>=clrsurf) { /* note the inversion! */
            lam_f=lv;
          } else {
            lam_f=ll;
            //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
          }
          real tmp_f;
          real len_f = distance_z(i,j,k,-1,tmp_f);
          len_f = 0.5*phi.dzc(k) - len_f;

          tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
        }
        assert(boil::realistic(tmp_w));
        
        tnl[i][j][k] = lam_s/ll*(tmp_w-tmp_s)/len_s * (-1);
      }
    }
  }
#endif

  for_ijk(i,j,k){
    if(Interface(i,j,k)){
      real lv = lambdav;
      real ll = lambdal;
      if (diff_eddy) {
        lv += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
        ll += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real qv = -lv*tnv[i][j][k];
      real ql =  ll*tnl[i][j][k];
      M[i][j][k] = (qv + ql) / latent;

#if 0
      if(i==3&&j==3&&k==3) {
      boil::aout<<"PCV_m: "<<i<<" "<<j<<" "<<k<<" | "<<phi.xc(i)<<" | "<<clr[i][j][k]<<" "<<adens[i][j][k]<<" "<<fs[Comp::k()][i][j][k]<<" "<<fs[Comp::k()][i+1][j][k]<<" | "<<tpr[i][j][k]<<" "<<tpr[i][j][k-1]<<" "<<tpr[i][j][k+1]<<" | "<<tpr[i-1][j][k]<<" "<<tpr[i+1][j][k]<<" | "<<qv<<" "<<ql<<" "<<qv+ql<<" | "<<tzv[i][j][k]<<" "<<nx[i][j][k]<<" "<<nz[i][j][k]<<boil::endl;
      }
#endif
    } else {
      M[i][j][k] = 0.0;
    }
  }

  M.bnd_update();
  M.exchange_all();

  return;
}

