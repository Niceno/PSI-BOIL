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
    inline void set_mass_resistance(const real mresisnew) {
      mresis = mresisnew;
      boil::oout<<"Schrage: Mass transfer resistance= "<<mresisnew<<boil::endl;
    }

    const Matter * fluid() const {return flu;}

  protected:
    real mresis, rhol;
    Matter * flu;

    const Scalar mflx;
    const Scalar * dpres;

    virtual void model();
    void pressure_effect();
    void mass_src_effect();
};

#endif
