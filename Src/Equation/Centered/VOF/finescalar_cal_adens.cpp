#include "finescalar.h"

/******************************************************************************/
void FineScalar::cal_adens() {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adens
*******************************************************************************/

  boil::timer.start("finescalar adens");

  /* cell centered */
  for_vijk(adens,i,j,k) {
    real gradient = 0.0;
   
    real dx = adens.dxc(i);
    real dy = adens.dyc(j);
    real dz = adens.dzc(k);
 
    /* x-gradient */
    real marker_w   = value(i,j,k, w  ()); 
    real marker_ws  = value(i,j,k, ws ()); 
    real marker_wn  = value(i,j,k, wn ()); 
    real marker_wb  = value(i,j,k, wb ()); 
    real marker_wt  = value(i,j,k, wt ()); 
    real marker_wsb = value(i,j,k, wsb()); 
    real marker_wst = value(i,j,k, wst()); 
    real marker_wnb = value(i,j,k, wnb()); 
    real marker_wnt = value(i,j,k, wnt()); 

    real marker_e   = value(i,j,k, e  ()); 
    real marker_es  = value(i,j,k, es ()); 
    real marker_en  = value(i,j,k, en ()); 
    real marker_eb  = value(i,j,k, eb ()); 
    real marker_et  = value(i,j,k, et ()); 
    real marker_esb = value(i,j,k, esb()); 
    real marker_est = value(i,j,k, est()); 
    real marker_enb = value(i,j,k, enb()); 
    real marker_ent = value(i,j,k, ent()); 

    real gradx = grad_1D(marker_w,
                         marker_ws, marker_wn, marker_wb, marker_wt,
                         marker_wsb, marker_wst, marker_wnb, marker_wnt,
                         marker_e,
                         marker_es, marker_en, marker_eb, marker_et,
                         marker_esb, marker_est, marker_enb, marker_ent,
                         dx, k);

    /* y-gradient */
    real marker_s   = value(i,j,k, s  ()); 
    real marker_sw  = marker_ws; 
    real marker_se  = marker_es; 
    real marker_sb  = value(i,j,k, sb ()); 
    real marker_st  = value(i,j,k, st ()); 
    real marker_swb = marker_wsb; 
    real marker_swt = marker_wst; 
    real marker_seb = marker_esb; 
    real marker_set = marker_est; 

    real marker_n   = value(i,j,k, n  ()); 
    real marker_nw  = marker_wn; 
    real marker_ne  = marker_en; 
    real marker_nb  = value(i,j,k, nb ()); 
    real marker_nt  = value(i,j,k, nt ()); 
    real marker_nwb = marker_wnb; 
    real marker_nwt = marker_wnt; 
    real marker_neb = marker_enb; 
    real marker_net = marker_ent; 

    real grady = grad_1D(marker_s,
                         marker_sw, marker_se, marker_sb, marker_st,
                         marker_swb, marker_swt, marker_seb, marker_set,
                         marker_n,
                         marker_nw, marker_ne, marker_nb, marker_nt,
                         marker_nwb, marker_nwt, marker_neb, marker_net,
                         dy, k);

    /* z-gradient */
    real marker_b   = value(i,j,k, b  ()); 
    real marker_bw  = marker_wb; 
    real marker_be  = marker_eb; 
    real marker_bs  = marker_sb; 
    real marker_bn  = marker_nb; 
    real marker_bws = marker_wsb; 
    real marker_bwn = marker_wnb; 
    real marker_bes = marker_esb; 
    real marker_ben = marker_enb; 

    real marker_t   = value(i,j,k, t  ()); 
    real marker_tw  = marker_wt; 
    real marker_te  = marker_et; 
    real marker_ts  = marker_st; 
    real marker_tn  = marker_nt; 
    real marker_tws = marker_wst; 
    real marker_twn = marker_wnt; 
    real marker_tes = marker_est; 
    real marker_ten = marker_ent; 

    real gradz = grad_1D(marker_b,
                         marker_bw, marker_be, marker_bs, marker_bn, 
                         marker_bws, marker_bwn, marker_bes, marker_ben,
                         marker_t,
                         marker_tw, marker_te, marker_ts, marker_tn,
                         marker_tws, marker_twn, marker_tes, marker_ten,
                         dz, k);

    adens[i][j][k] = sqrt(gradx*gradx+grady*grady+gradz*gradz);
    //adens[i][j][k] = ((*phi)[i][j][k]>0.)*sqrt(gradx*gradx+grady*grady+gradz*gradz);
    //adens[i][j][k]=2.*(*phi)[i][j][k]*sqrt(gradx*gradx+grady*grady+gradz*gradz);
    //adens[i][j][k]=6.*(1.-(*phi)[i][j][k])*(*phi)[i][j][k]*sqrt(gradx*gradx+grady*grady+gradz*gradz);

  }

  real sum(0.0);
  int count(0);
  for_vijk(adens,i,j,k) {
    real sumplus =  adens[i][j][k]*adens.dV(i,j,k);
    if(sumplus>boil::atto) count++;
    sum += sumplus;
  }
  boil::oout<<"VOF::finescalar_adens "<<count<<" "<<sum<<boil::endl;

  boil::timer.stop("finescalar adens");

  return;
}

/*---------------------+
|  ancillary function  |
+---------------------*/
real FineScalar::grad_1D(const real mF,
                         const real mE1, const real mE2,
                         const real mE3, const real mE4,
                         const real mN1, const real mN2,
                         const real mN3, const real mN4,
                         const real pF,
                         const real pE1, const real pE2,
                         const real pE3, const real pE4,
                         const real pN1, const real pN2,
                         const real pN3, const real pN4,
                         const real delta, const int k) {
  /* no division by zero */
  if(delta<boil::atto) return 0.0;

  real wF = 4.0;
  real wE = 2.0;
  real wN = 1.0;

  real wT = wF*1.0 + wE*4.0 + wN*4.0;

  /* face */
  real gradval = wF*(pF-mF);
  
  /* edges */
  gradval += wE*(pE4-mE4 + pE3-mE3 + pE2-mE2 + pE1-mE1);

  /* nodes */
  gradval += wN*(pN4-mN4 + pN3-mN3 + pN2-mN2 + pN1-mN1);

  gradval /= wT;

  gradval /= delta;

#if 0
  if(fabs(pF-mF)>0.0) {
    boil::oout<<"HERE! "<<k<<" : "<<pF<<" "<<mF<<" "<<pE4<<" "<<mE4<<" "<<pE3<<" "<<mE3<<" "<<pE2<<" "<<mE2<<" "<<pE1<<" "<<mE1<<" "<<pN4<<" "<<mN4<<" "<<pN3<<" "<<mN3<<" "<<pN2<<" "<<mN2<<" "<<pN1<<" "<<mN1<<boil::endl;
  }
#endif

  return gradval;
}

