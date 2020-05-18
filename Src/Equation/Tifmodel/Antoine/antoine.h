#ifndef ANTOINE_H
#define ANTOINE_H

#include "../tif.h"


///////////////////////////////////////////
//                                       //
//  Antoine interface temperature model  //
//                                       //
///////////////////////////////////////////
/* constants should be supplied for pressure in mmHg, as standard */
class Antoine : public TIF {
  public:
    Antoine(const real tref,
            const Topology * topo,
            const Scalar & eps,
            const real A,
            const real B,
            const real C);
    ~Antoine() {}

    real value(const int i, const int j, const int k);
    real temperature(const real eps);
    real epsilon(const real tint);
    real pressure(const real tint);

  protected:
    const Scalar eps;
    real A,B,C,C_K;

    real an_value(const real alpha);

    virtual void model();
};

#endif
