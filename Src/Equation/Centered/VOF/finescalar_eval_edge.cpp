#include "finescalar.h"

/******************************************************************************/
void FineScalar::eval_edge() {
/***************************************************************************//**
*  \brief Evaluate the marker function at edges.
*******************************************************************************/

  /* staggered in y and z */
  for(int i=phi->si(); i<=phi->ei()  ; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i][j-1][k-1];
    real phi_mp = (*phi)[i][j-1][k  ];
    real phi_pm = (*phi)[i][j  ][k-1];
    real phi_pp = (*phi)[i][j  ][k  ];

    /* evaluate each cell pair */
#if 0
    bool mm_mp = (phi_mm-phisurf)*(phi_mp-phisurf)>0.0;
    bool mm_pm = (phi_mm-phisurf)*(phi_pm-phisurf)>0.0;
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;
    bool mp_pp = (phi_mp-phisurf)*(phi_pp-phisurf)>0.0;
    bool pm_pp = (phi_pm-phisurf)*(phi_pp-phisurf)>0.0;
#else
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;
#endif

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mm, i,j-1,k-1, nt(),
                              phi_pp, i,j  ,k  , sb());
    }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mp, i,j-1,k  , nb(),
                              phi_pm, i,j  ,k-1, st());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,sb()) = sum/2.0;
     
  }

  /* staggered in x and z */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()  ; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i-1][j][k-1];
    real phi_mp = (*phi)[i-1][j][k  ];
    real phi_pm = (*phi)[i  ][j][k-1];
    real phi_pp = (*phi)[i  ][j][k  ];

    /* evaluate each cell pair */
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mm, i-1,j,k-1, et(),
                              phi_pp, i  ,j,k  , wb());
    }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mp, i-1,j,k  , eb(),
                              phi_pm, i  ,j,k-1, wt());
    }

    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,wb()) = sum/2.0;
     
  }

  /* staggered in x and y */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()  ; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i-1][j-1][k];
    real phi_mp = (*phi)[i-1][j  ][k];
    real phi_pm = (*phi)[i  ][j-1][k];
    real phi_pp = (*phi)[i  ][j  ][k];

    /* evaluate each cell pair */
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mm, i-1,j-1,k, en(),
                              phi_pp, i  ,j  ,k, ws());
     }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mp, i-1,j  ,k, es(),
                              phi_pm, i  ,j-1,k, wn());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,ws()) = sum/2.0;
     
  }

  /* correct at boundaries */
  //bdcond_edge();

  return;
}
