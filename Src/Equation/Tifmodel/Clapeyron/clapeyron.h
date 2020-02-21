#ifndef CLAPEYRON_H
#define CLAPEYRON_H

#include "../tif.h"


//////////////////////////////////////////////////////
//                                                  //
//  Clausius-Clapeyron interface temperature model  //
//                                                  //
//////////////////////////////////////////////////////
class Clapeyron : public TIF {
  public:
    Clapeyron(const real tref,
              const Scalar * adens,
              const Scalar & eps,
              const real mv,
              const real latent,
              const real latent_slp = 0.0);
    ~Clapeyron() {}

    real value(const int i, const int j, const int k);
    real temperature(const real eps);
    real epsilon(const real tint);

    inline int get_iteration_count() const { return nmax; }
    inline real get_iteration_tol() const  { return errmax; }
    inline void set_iteration_params(const int nnew, const real enew) {
      nmax = nnew;
      errmax = enew;
      boil::oout<<"Clapeyron::iterationparams: "<<nnew<<" "<<enew<<"\n";
    }

  protected:
    const Scalar eps;
    real tri;
    real latent, latent_cst, latent_slp, Rm;

    real errmax;
    int nmax;

    virtual void model();

    real cc_iterate(const real alpha);
    real err(const real tau, const real alpha);
    real err_prime(const real tau);
};

#endif
