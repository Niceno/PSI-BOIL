#ifndef SITE_H
#define SITE_H

#include <vector>
#include "../../Global/global_precision.h"
#include "../../Global/global_constants.h"
#include "../../Global/global_realistic.h"

////////////
//        //
//  Site  //
//        //
////////////
class Site {

  public:
    Site(const real x, const real y, const real z
       , const real t1, const real t2);
    ~Site();

    inline real x() const {return xsite;};
    inline real y() const {return ysite;};
    inline real z() const {return zsite;};
    inline int ic_site() const {return icsite;}
    inline int jc_site() const {return jcsite;}
    inline int kc_site() const {return kcsite;}

    inline int is() const {return ist;};
    inline int ie() const {return ied;};
    inline int js() const {return jst;};
    inline int je() const {return jed;};
    inline int ks() const {return kst;};
    inline int ke() const {return ked;};

    void set_x (real r){xsite=r;};
    void set_y (real r){ysite=r;};
    void set_z (real r){zsite=r;};

    void set_ic_site (int i){icsite=i;};
    void set_jc_site (int j){jcsite=j;};
    void set_kc_site (int k){kcsite=k;};

    void set_is(int i) {ist=i;};
    void set_ie(int i) {ied=i;};
    void set_js(int j) {jst=j;};
    void set_je(int j) {jed=j;};
    void set_ks(int k) {kst=k;};
    void set_ke(int k) {ked=k;};

    /* time_Tact: time when Tw > Tact */
    void set_time_Tact(real t) {tm_Tact=t;};
    real time_Tact() {return tm_Tact;}

    /* time_plant_clr: time when color function for bubble is generated */
    void set_time_plant_clr(real t) {tm_plant_clr=t;};
    real time_plant_clr() {return tm_plant_clr;}

    /* time_cutneck */
    void set_time_cutneck(real t) {tcutneck=t;};
    real time_cutneck() {return tcutneck;}

    /* plant color function in current time step */
    void set_plant_clr_current(bool b) {b_plant_clr_current=b;};
    bool plant_clr_current() {return b_plant_clr_current;}

    /* plant color function in current time step */
    void set_heat_sink_current(bool b) {b_heat_sink_current=b;};
    bool heat_sink_current() {return b_heat_sink_current;}

    /* neck */
    void set_neck(bool b) {bneck=b;};
    bool neck() {return bneck;}

    /* neck previous time step */
    void set_neck_prev(bool b) {bneck_prev=b;};
    bool neck_prev() {return bneck_prev;}

    /* alow_replant from outside */
    void set_allow_replant(bool b) {allow_seed=b;};
    bool allow_replant() {return allow_seed;};

    /* request to outside to allow replant */
    void set_req_replant(bool b) {req_seed=b;};
    bool req_replant() {return req_seed;};

    /* active temperature */
    void set_active_tpr(real r) {act_tpr=r;};
    real active_tpr() {return act_tpr;};

    /* height of replant */
    void set_zplant(real r) {zplt=r;};
    real zplant() {return zplt;};

    /* fater: parent ID of dummy site */
    void set_father(int i) {fthr=i;};
    int father() {return fthr;};

    /* active: nucleation site is active 
     * from the time Tw > Tact (plant/replant)
     * till the time, when color at nucleation site is covered with liquid */
    void set_active(const bool b) {act=b;}
    bool active() const {return act;}

    /* qsink: add heat sink to qsrc */
    void set_qsink(const bool b) {bqsink=b;}
    bool qsink() {return bqsink;}

    /* bubble volume */
    void set_vol_bubble(real r) {v_bubble=r;}
    real vol_bubble() {return v_bubble;}

    /* bubble base area */
    void set_area_base(real r) {a_base=r;}
    real area_base() {return a_base;}

    /* sum sink energy */
    void set_sink_energy(real r) {s_energy=r;}
    real sink_energy() {return s_energy;}

    /* sum sink energy in time */
    void set_sum_sink_energy(real r) {s_sum_energy=r;}
    real sum_sink_energy() {return s_sum_energy;}

    /* k-adjacent in fluid side */
    void set_kadj(int i) {kadjacent=i;}
    real kadj() {return kadjacent;}

    /* contain center of bubble
     * true: center of bubble is inside of decomposed domain
     * false: center of bubble is outside of decomposed domain */
    //void set_contain_cb(bool b) {containCb=b;}
    //bool contain_cb() {return containCb;}

    /* contain nucleation site (z=0)
     * true: site is inside of decomposed domain
     * false: site is outside of decomposed domain */
    void set_contain_site(bool b) {containSite=b;}
    bool contain_site() {return containSite;}

    /* contain range (loop range around site)
     * true: range is inside of decomposed domain
     * false: range is outside of decomposed domain */
    void set_contain_range(bool b) {containRange=b;}
    bool contain_range() {return containRange;}

    /* UNUSED: seed previous time step */
    void set_seed_prev(bool b) {bseed_prev=b;};
    bool seed_prev() {return bseed_prev;}

  private:
    real xsite, ysite, zsite;
    real tm_Tact, tm_plant_clr, tcutneck;
    real act_tpr, zplt;
    int ist, ied, jst, jed, kst, ked;
    int icsite, jcsite, kcsite, kadjacent;
    bool b_plant_clr_current, b_heat_sink_current;
    bool bneck, bneck_prev;
    bool req_seed, allow_seed;
    bool act;
    bool bqsink, containCb, containRange, containSite;
    int fthr;
    real v_bubble, a_base;
    real s_energy,s_sum_energy;
    // UNUSED
    bool bseed_prev;
};

#endif
