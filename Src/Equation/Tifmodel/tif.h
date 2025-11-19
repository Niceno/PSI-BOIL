#ifndef TIF_H
#define TIF_H

#include "../../Field/Scalar/scalar.h"
#include "../../Matter/matter.h"
#include "../../Global/global_realistic.h"

///////////////////////////////////
//                               //
//  Interface temperature model  //
//                               //
///////////////////////////////////
class TIF {
  public:
    TIF(const real tref); 
    TIF(const real tref, 
        const real latent,
        const real mresis,
        Matter * flu,
        const Scalar & adens,
        const Scalar & mflx,
        const Scalar * pres = NULL);
    ~TIF() {}

    void update_tifold() {
      for_vijk(tifold,i,j,k) 
        tifold[i][j][k] = tif[i][j][k];
      store_tif = true;
      tifold.exchange();
    }
    void tint_field(const bool newstep = true);

    real get_ur() { return factor; }
    void set_ur(const real factnew) {
      factor = factnew;
      boil::oout<<"TIFmodel: Underrelaxation factor = "<<factnew<<boil::endl;
    }

    real tref() { return tr; }
    void set_tref(const real tnew) {
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
    void set_strong_limiting(const Scalar * tpr,
                             const Scalar * clr,
                             const real clrsurf);

  protected:
    real factor; /* under-relaxation factor */
    real tmin, tmax;
    bool weaklim, stronglim;
    real clrsurf;

    bool store_tif,variable_tif;
    real tr, latent, mresis, rhol;

    const Matter * fluid() const {return flu;}
    Matter * flu;

    const Scalar mflx;
    const Scalar adens;
    const Scalar * dpres;
    const Scalar * clr;
    const Scalar * tpr;

    void Pressure_effect();
    void Mass_src_effect();
    void Extend_tint    ();

    bool Vicinity(const int i, const int j, const int k);
    bool Interface(const int i, const int j, const int k);
    inline real Underrelaxation(const real tintnew, const real tifold);
};

#endif
