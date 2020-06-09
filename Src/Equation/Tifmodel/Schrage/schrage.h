#ifndef SCHRAGE_H
#define SCHRAGE_H

#include "../tif.h"

///////////////////////////////////////////
//                                       //
//  Schrage interface temperature model  //
//      (with disjoining pressure)       //
//                                       //
///////////////////////////////////////////
class Schrage : public TIF {
  public:
    Schrage(const real tref, 
            Matter * flu,
            Topology * topo,
            const Scalar & mflx,
            const Scalar * pres = NULL);
    ~Schrage() {}

    inline real get_mass_resistance() { return mresis; }
    inline real get_heat_resistance() { return hresis; }
    inline void set_heat_resistance(const real resisnew) {
      hresis = resisnew;
      mresis = hresis*latent;
      boil::oout<<"Schrage: Transfer resistance= "<<hresis
                <<" "<<mresis<<boil::endl;
    }

    const Matter * fluid() const {return flu;}

    static real calculate_heat_transfer_resistance(
                                          const real tr, const real rhov,
                                          const real mmass, const real latent);

  protected:
    real mresis, hresis, latent, rhol;
    Matter * flu;

    const Scalar mflx;
    const Scalar * dpres;

    virtual void model();
    void pressure_effect();
    void mass_src_effect();
};

#endif
