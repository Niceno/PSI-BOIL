#include "finescalar.h"

/******************************************************************************/
void FineScalar::eval_node() {
/***************************************************************************//**
*  \brief Evaluate the marker function at nodes.
*******************************************************************************/

  /* staggered in x and y and z */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mmm = (*phi)[i-1][j-1][k-1];
    real phi_mmp = (*phi)[i-1][j-1][k  ];
    real phi_mpm = (*phi)[i-1][j  ][k-1];
    real phi_mpp = (*phi)[i-1][j  ][k  ];
    real phi_pmm = (*phi)[i  ][j-1][k-1];
    real phi_pmp = (*phi)[i  ][j-1][k  ];
    real phi_ppm = (*phi)[i  ][j  ][k-1];
    real phi_ppp = (*phi)[i  ][j  ][k  ];

    /* evaluate each cell pair */
    bool mmm_ppp = (phi_mmm-phisurf)*(phi_ppp-phisurf)>0.0;
    bool mmp_ppm = (phi_mmp-phisurf)*(phi_ppm-phisurf)>0.0;
    bool mpp_pmm = (phi_mpp-phisurf)*(phi_pmm-phisurf)>0.0;
    bool mpm_pmp = (phi_mpm-phisurf)*(phi_pmp-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[4];
    int idx(0);

    if(mmm_ppp) {
      vals[idx] = phi_mmm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mmm, i-1,j-1,k-1, ent(),
                              phi_ppp, i  ,j  ,k  , wsb());
    }
    idx++;
    if(mmp_ppm) {
      vals[idx] = phi_mmp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mmp, i-1,j-1,k  , enb(),
                              phi_ppm, i  ,j  ,k-1, wst());
    }
    idx++;
    if(mpp_pmm) {
      vals[idx] = phi_mpp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mpp, i-1,j  ,k  , esb(),
                              phi_pmm, i  ,j-1,k-1, wnt());
    }
    idx++;
    if(mpm_pmp) {
      vals[idx] = phi_mpm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mpm, i-1,j  ,k-1, est(),
                              phi_pmp, i  ,j-1,k  , wnb());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k, wsb()) = sum/4.0;
  }

  /* correct at boundaries */
  //bdcond_node();

  return;
}
