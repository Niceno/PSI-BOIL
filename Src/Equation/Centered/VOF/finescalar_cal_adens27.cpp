#include "finescalar.h"

/*----------------------+
|  ancillary functions  |
+----------------------*/
void cal_1D_derivatives(real & grad_00_m, real & grad_00_p,
                        real & grad_m0_m, real & grad_m0_p,
                        real & grad_p0_m, real & grad_p0_p,
                        real & grad_0m_m, real & grad_0m_p,
                        real & grad_0p_m, real & grad_0p_p,
                        real & grad_mm_m, real & grad_mm_p,
                        real & grad_mp_m, real & grad_mp_p,
                        real & grad_pm_m, real & grad_pm_p,
                        real & grad_pp_m, real & grad_pp_p,
                        const real marker_000, const real marker_m00, const real marker_p00,
                        const real marker_0m0, const real marker_mm0, const real marker_pm0,
                        const real marker_0p0, const real marker_mp0, const real marker_pp0,
                        const real marker_00m, const real marker_m0m, const real marker_p0m,
                        const real marker_00p, const real marker_m0p, const real marker_p0p,
                        const real marker_0mm, const real marker_mmm, const real marker_pmm,
                        const real marker_0mp, const real marker_mmp, const real marker_pmp,
                        const real marker_0pm, const real marker_mpm, const real marker_ppm,
                        const real marker_0pp, const real marker_mpp, const real marker_ppp,
                        const real delta) {

    grad_00_m = (marker_000 - marker_m00)/delta;
    grad_00_p = (marker_p00 - marker_000)/delta;

    grad_m0_m = (marker_0m0 - marker_mm0)/delta;
    grad_m0_p = (marker_pm0 - marker_0m0)/delta;

    grad_p0_m = (marker_0p0 - marker_mp0)/delta;
    grad_p0_p = (marker_pp0 - marker_0p0)/delta;

    grad_0m_m = (marker_00m - marker_m0m)/delta;
    grad_0m_p = (marker_p0m - marker_00m)/delta;

    grad_0p_m = (marker_00p - marker_m0p)/delta;
    grad_0p_p = (marker_p0p - marker_00p)/delta;

    grad_mm_m = (marker_0mm - marker_mmm)/delta;
    grad_mm_p = (marker_pmm - marker_0mm)/delta;

    grad_mp_m = (marker_0mp - marker_mmp)/delta;
    grad_mp_p = (marker_pmp - marker_0mp)/delta;

    grad_pm_m = (marker_0pm - marker_mpm)/delta;
    grad_pm_p = (marker_ppm - marker_0pm)/delta;

    grad_pp_m = (marker_0pp - marker_mpp)/delta;
    grad_pp_p = (marker_ppp - marker_0pp)/delta;

    return;
}

inline real average4(const real a1, const real a2, const real a3, const real a4) {
    return 0.25*(a1+a2+a3+a4);
}

void cal_oct_derivatives(real & grad_mmm, real & grad_pmm, 
                         real & grad_mmp, real & grad_pmp,
                         real & grad_mpm, real & grad_ppm,
                         real & grad_mpp, real & grad_ppp,
                         const real grad_00_m, const real grad_00_p,
                         const real grad_m0_m, const real grad_m0_p,
                         const real grad_p0_m, const real grad_p0_p,
                         const real grad_0m_m, const real grad_0m_p,
                         const real grad_0p_m, const real grad_0p_p,
                         const real grad_mm_m, const real grad_mm_p,
                         const real grad_mp_m, const real grad_mp_p,
                         const real grad_pm_m, const real grad_pm_p,
                         const real grad_pp_m, const real grad_pp_p) {
    
    grad_mmm = average4(grad_mm_m, grad_m0_m, grad_0m_m, grad_00_m);
    grad_pmm = average4(grad_mm_p, grad_m0_p, grad_0m_p, grad_00_p);
    grad_mmp = average4(grad_mp_m, grad_m0_m, grad_0p_m, grad_00_m);
    grad_pmp = average4(grad_mp_p, grad_m0_p, grad_0p_p, grad_00_p);
    grad_mpm = average4(grad_pm_m, grad_p0_m, grad_0m_m, grad_00_m);
    grad_ppm = average4(grad_pm_p, grad_p0_p, grad_0m_p, grad_00_p);
    grad_mpp = average4(grad_pp_m, grad_p0_m, grad_0p_m, grad_00_m);
    grad_ppp = average4(grad_pp_p, grad_p0_p, grad_0p_p, grad_00_p);

    return;
}

inline real adensval(const real a1, const real a2, const real a3) {
    return sqrt(a1*a1+a2*a2+a3*a3);
}

/******************************************************************************/
void FineScalar::cal_adens27() {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adens
*******************************************************************************/

  boil::timer.start("finescalar adens27");

  /* cell centered */
  for_vijk(adens27,i,j,k) {
    real gradient = 0.0;
   
    real dx = adens27.dxc(i)/2.0;
    real dy = adens27.dyc(j)/2.0;
    real dz = adens27.dzc(k)/2.0;
 
    real marker_c;
    real color_c = (*phi)[i][j][k];
  
    if(color_c>phisurf)      marker_c = 1.0;
    else if(color_c<phisurf) marker_c = 0.0;
    else                     marker_c = 0.5;

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

    /* x-gradients */
    real xgrad_c_m, xgrad_c_p, xgrad_s_m, xgrad_s_p, xgrad_n_m, xgrad_n_p,
         xgrad_b_m, xgrad_b_p, xgrad_t_m, xgrad_t_p, xgrad_sb_m, xgrad_sb_p,
         xgrad_st_m, xgrad_st_p, xgrad_nb_m, xgrad_nb_p, xgrad_nt_m, xgrad_nt_p;

    cal_1D_derivatives(xgrad_c_m, xgrad_c_p, xgrad_s_m, xgrad_s_p,
                       xgrad_n_m, xgrad_n_p, xgrad_b_m, xgrad_b_p,
                       xgrad_t_m, xgrad_t_p, xgrad_sb_m, xgrad_sb_p,
                       xgrad_st_m, xgrad_st_p, xgrad_nb_m, xgrad_nb_p,
                       xgrad_nt_m, xgrad_nt_p,
                       marker_c, marker_w, marker_e, 
                       marker_s, marker_ws, marker_es,
                       marker_n, marker_wn, marker_en,
                       marker_b, marker_wb, marker_eb,
                       marker_t, marker_wt, marker_et,
                       marker_sb, marker_wsb, marker_esb,
                       marker_st, marker_wst, marker_est,
                       marker_nb, marker_wnb, marker_enb,
                       marker_nt, marker_wnt, marker_ent,
                       dx);

    /* y-gradients */
    real ygrad_c_m, ygrad_c_p, ygrad_w_m, ygrad_w_p, ygrad_e_m, ygrad_e_p,
         ygrad_b_m, ygrad_b_p, ygrad_t_m, ygrad_t_p, ygrad_wb_m, ygrad_wb_p,
         ygrad_wt_m, ygrad_wt_p, ygrad_eb_m, ygrad_eb_p, ygrad_et_m, ygrad_et_p;

   cal_1D_derivatives(ygrad_c_m, ygrad_c_p, ygrad_w_m, ygrad_w_p,
                      ygrad_e_m, ygrad_e_p, ygrad_b_m, ygrad_b_p,
                      ygrad_t_m, ygrad_t_p,  ygrad_wb_m, ygrad_wb_p,
                      ygrad_wt_m, ygrad_wt_p,  ygrad_eb_m, ygrad_eb_p,
                      ygrad_et_m, ygrad_et_p,
                      marker_c, marker_s, marker_n,
                      marker_w, marker_sw, marker_nw,
                      marker_e, marker_se, marker_ne,
                      marker_b, marker_sb, marker_nb,
                      marker_t, marker_st, marker_nt,
                      marker_wb, marker_swb, marker_nwb,
                      marker_wt, marker_swt, marker_nwt,
                      marker_eb, marker_seb, marker_neb,
                      marker_et, marker_set, marker_net,
                      dy);

    /* z-gradients */
    real zgrad_c_m, zgrad_c_p, zgrad_w_m, zgrad_w_p, zgrad_e_m, zgrad_e_p,
         zgrad_s_m, zgrad_s_p, zgrad_n_m, zgrad_n_p, zgrad_ws_m, zgrad_ws_p, 
         zgrad_es_m, zgrad_es_p, zgrad_wn_m, zgrad_wn_p, zgrad_en_m, zgrad_en_p;

    cal_1D_derivatives(zgrad_c_m, zgrad_c_p, zgrad_w_m, zgrad_w_p,
                       zgrad_e_m, zgrad_e_p, zgrad_s_m, zgrad_s_p,
                       zgrad_n_m, zgrad_n_p, zgrad_ws_m, zgrad_ws_p,
                       zgrad_wn_m, zgrad_wn_p, zgrad_es_m, zgrad_es_p,
                       zgrad_en_m, zgrad_en_p,
                       marker_c, marker_b, marker_t,
                       marker_w, marker_bw, marker_tw,
                       marker_e, marker_be, marker_te,
                       marker_s, marker_bs, marker_ts,
                       marker_n, marker_bn, marker_tn,
                       marker_ws, marker_bws, marker_tws,
                       marker_wn, marker_bwn, marker_twn,
                       marker_es, marker_bes, marker_tes,
                       marker_en, marker_ben, marker_ten,
                       dz);
    
    /* x octants */
    real xgrad_wsb, xgrad_esb, xgrad_wst, xgrad_est,
         xgrad_wnb, xgrad_enb, xgrad_wnt, xgrad_ent;

    cal_oct_derivatives(xgrad_wsb, xgrad_esb, xgrad_wst, xgrad_est,
                        xgrad_wnb, xgrad_enb, xgrad_wnt, xgrad_ent,
                        xgrad_c_m, xgrad_c_p, xgrad_s_m, xgrad_s_p,
                        xgrad_n_m, xgrad_n_p, xgrad_b_m, xgrad_b_p,
                        xgrad_t_m, xgrad_t_p, xgrad_sb_m, xgrad_sb_p,
                        xgrad_st_m, xgrad_st_p, xgrad_nb_m, xgrad_nb_p,
                        xgrad_nt_m, xgrad_nt_p);

    /* y octants */
    real ygrad_swb, ygrad_nwb, ygrad_swt, ygrad_nwt,
         ygrad_seb, ygrad_neb, ygrad_set, ygrad_net;

    cal_oct_derivatives(ygrad_swb, ygrad_nwb, ygrad_swt, ygrad_nwt,
                        ygrad_seb, ygrad_neb, ygrad_set, ygrad_net,
                        ygrad_c_m, ygrad_c_p, ygrad_w_m, ygrad_w_p,
                        ygrad_e_m, ygrad_e_p, ygrad_b_m, ygrad_b_p,
                        ygrad_t_m, ygrad_t_p, ygrad_wb_m, ygrad_wb_p,
                        ygrad_wt_m, ygrad_wt_p, ygrad_eb_m, ygrad_eb_p,
                        ygrad_et_m, ygrad_et_p);

    /* z octants */
    real zgrad_bws, zgrad_tws, zgrad_bwn, zgrad_twn,
         zgrad_bes, zgrad_tes, zgrad_ben, zgrad_ten;

    cal_oct_derivatives(zgrad_bws, zgrad_tws, zgrad_bwn, zgrad_twn,
                        zgrad_bes, zgrad_tes, zgrad_ben, zgrad_ten,
                        zgrad_c_m, zgrad_c_p, zgrad_w_m, zgrad_w_p,
                        zgrad_e_m, zgrad_e_p, zgrad_s_m, zgrad_s_p,
                        zgrad_n_m, zgrad_n_p, zgrad_ws_m, zgrad_ws_p,
                        zgrad_wn_m, zgrad_wn_p, zgrad_es_m, zgrad_es_p,
                        zgrad_en_m, zgrad_en_p);
     
    /* xyz octants */
    real area_wsb, area_esb, area_wst, area_est,
         area_wnb, area_enb, area_wnt, area_ent;

    real vol   = (*phi).dV(i,j,k);
    real voloct= vol/8.0;
  
    area_wsb = adensval(xgrad_wsb,ygrad_swb,zgrad_bws);
    area_esb = adensval(xgrad_esb,ygrad_seb,zgrad_bes);
    area_wnb = adensval(xgrad_wnb,ygrad_nwb,zgrad_bwn);
    area_enb = adensval(xgrad_enb,ygrad_neb,zgrad_ben);
    area_wst = adensval(xgrad_wst,ygrad_swt,zgrad_tws);
    area_est = adensval(xgrad_est,ygrad_set,zgrad_tes);
    area_wnt = adensval(xgrad_wnt,ygrad_nwt,zgrad_twn);
    area_ent = adensval(xgrad_ent,ygrad_net,zgrad_ten);

    area_wsb *= voloct;
    area_esb *= voloct;
    area_wnb *= voloct;
    area_enb *= voloct;
    area_wst *= voloct;
    area_est *= voloct;
    area_wnt *= voloct;
    area_ent *= voloct;

    /* cell */
    real area_c = area_wsb+area_esb+area_wnb+area_enb
                 +area_wst+area_est+area_wnt+area_ent;


    adens27[i][j][k] = area_c/vol;

  }

  real sum(0.0);
  int count(0);
  for_vijk(adens27,i,j,k) {
    real sumplus =  adens27[i][j][k]*adens27.dV(i,j,k);
    if(sumplus>boil::atto) count++;
    sum += sumplus;
  }
  boil::oout<<"VOF::finescalar_adens27 "<<count<<" "<<sum<<boil::endl;

  boil::timer.stop("finescalar adens27");

  return;
}
