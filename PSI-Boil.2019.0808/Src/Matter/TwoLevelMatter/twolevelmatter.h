#ifndef TWOLEVELMATTER_H
#define TWOLEVELMATTER_H

#include "../matter.h"

////////////////////////
//                    //
//  Two-level Matter  //
//                    //
////////////////////////
class TwoLevelMatter {
  public:
    TwoLevelMatter(const TwoLevelDomain & d,
                   const char * nm = NULL) :
      coarse(d.coarse(),nm),
      fine(d.fine(),nm)
    {
    }

    TwoLevelMatter(const TwoLevelMatter & a,
                   const TwoLevelMatter & b,
                   const TwoLevelScalar * ca) :       
      coarse(a.coarse,b.coarse,&ca->coarse),
      fine(a.fine,b.fine,&ca->fine)
    {
    }

    TwoLevelMatter(const TwoLevelMatter & a,
                   const TwoLevelMatter & b,
                   const TwoLevelScalar * ca,      
                   const TwoLevelVector * bdca) :     
      coarse(a.coarse,b.coarse,&ca->coarse,&bdca->coarse),
      fine(a.fine,b.fine,&ca->fine,&bdca->fine)
    {
    }

    /* set (initialize) physical properties */
    void rho   (const real & v) { coarse.rho(v);    fine.rho(v); }
    void mu    (const real & v) { coarse.mu(v);     fine.mu(v); }
    void cp    (const real & v) { coarse.cp(v);     fine.cp(v); }
    void lambda(const real & v) { coarse.lambda(v); fine.lambda(v); }
    void gamma (const real & v) { coarse.gamma(v);  fine.gamma(v); }
    void beta  (const real & v) { coarse.beta(v);   fine.beta(v); }
    void mmass (const real & v) { coarse.mmass(v);  fine.mmass(v); }
    void sigma (const real & v) { coarse.sigma(v);  fine.sigma(v); }
    void latent(const real & v) { coarse.latent(v); fine.latent(v); }

    /* rescale properties by factors of length and time */
    void rescale(const real xmult = 1., const real tmult = 1.,
                 const real mmult = 1.) {
       coarse.rescale(xmult,tmult,mmult);
       fine.rescale(xmult,tmult,mmult);
    }

    /* set certain property to be variable */
    void variable(const Set & s) {
       coarse.variable(s);
       fine.variable(s);
    }

    /* set all properties to be variable */
    void variable() {
       coarse.variable();
       fine.variable();
    }

    /* check if matter is mixture */
    inline bool mixture() const { return coarse.mixture(); }

    Matter coarse;
    Matter fine;
};
#endif
