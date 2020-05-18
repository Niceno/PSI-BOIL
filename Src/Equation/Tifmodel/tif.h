#ifndef TIF_H
#define TIF_H

#include "../../Field/Scalar/scalar.h"
#include "../../Field/ScalarInt/scalarint.h"
#include "../Topology/topology.h"
#include "../../Matter/matter.h"
#include "../../Global/global_realistic.h"
#include "../../Plot/plot.h"

///////////////////////////////////
//                               //
//  Interface temperature model  //
//                               //
///////////////////////////////////
class TIF {
  public:
    TIF(const real tref); 
    TIF(const real tref, Topology * topo); 
    ~TIF() {}

    void init() {
      if(variable_tif) { 
        model(); 
        tif.bnd_update();
        tif.exchange_all();
        extend_tint();
      }
    }

    void update_tifold() {
      tifold = tif;
      store_tif = true;
    }
    void tint_field(const bool newstep = true);

    inline real get_ur() const { return factor; }
    inline void set_ur(const real factnew) {
      factor = factnew;
      boil::oout<<"TIFmodel: Underrelaxation factor = "<<factnew<<boil::endl;
    }

    inline real tref() const { return tr; }
    inline void set_tref(const real tnew) {
      tr = tnew;
      boil::oout<<"TIFmodel: Tref = "<<tnew<<boil::endl;
    }

    Scalar tif, tifold;

    real Tint(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k) const;
    real Tint_old(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k) const;
    real Tint(const int i, const int j, const int k) const;
    real Tint_old(const int i, const int j, const int k) const; 

    void set_weak_limiting(const real tmin, const real tmax);

  protected:
    real factor; /* under-relaxation factor */
    real tmin, tmax;
    bool weaklim;

    bool store_tif,variable_tif;
    real tr;

    Topology * topo;
    ScalarInt iflag, tempflag, tempflag2;
    Scalar stmp;

    virtual void model() {};
    void extend_tint();

    inline real underrelaxation(const real tintnew, const real tifold);
};

#endif
