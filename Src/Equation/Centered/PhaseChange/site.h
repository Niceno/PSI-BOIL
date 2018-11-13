#ifndef SITE_H
#define SITE_H

#include <vector>
#include "../../../Global/global_precision.h"
#include "../../../Global/global_constants.h"

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

    real x() {return xsite;};
    real y() {return ysite;};
    real z() {return zsite;};
    int ic() {return icent;}
    int jc() {return jcent;}
    int kc() {return kcent;}

    int is() {return ist;};
    int ie() {return ied;};
    int js() {return jst;};
    int je() {return jed;};
    int ks() {return kst;};
    int ke() {return ked;};

    void set_x (real r){xsite=r;};
    void set_y (real r){ysite=r;};
    void set_z (real r){zsite=r;};

    void set_ic (int i){icent=i;};
    void set_jc (int j){jcent=j;};
    void set_kc (int k){kcent=k;};

    void set_is(int i) {ist=i;};
    void set_ie(int i) {ied=i;};
    void set_js(int j) {jst=j;};
    void set_je(int j) {jed=j;};
    void set_ks(int k) {kst=k;};
    void set_ke(int k) {ked=k;};

    /* time_seed */
    void set_time_seed(real t) {tseed=t;};
    real time_seed() {return tseed;}

    /* time_seed */
    void set_time_cutneck(real t) {tcutneck=t;};
    real time_cutneck() {return tcutneck;}

    /* seed */
    void set_seed(bool b) {bseed=b;};
    bool seed() {return bseed;}

    /* seed previous time step*/
    void set_seed_prev(bool b) {bseed_prev=b;};
    bool seed_prev() {return bseed_prev;}

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

    /* active: nucleation site is active or not */
    void set_active(const bool b) {act=b;}
    bool active() {return act;}

    /* qsink: add heat sink to qsrc */
    void set_qsink(const bool b) {bqsink=b;}
    bool qsink() {return bqsink;}

  private:
    real xsite, ysite, zsite;
    real tseed, tcutneck;
    real act_tpr, zplt;
    int ist, ied, jst, jed, kst, ked;
    int icent, jcent, kcent;
    bool bseed, bneck, bseed_prev, bneck_prev;
    bool req_seed, allow_seed;
    bool act;
    bool bqsink;
    int fthr;

};

#endif

