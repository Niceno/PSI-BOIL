#ifndef PROPERTY_H
#define PROPERTY_H

#include <vector>

#include "../Domain/domain.h"
#include "../Field/Scalar/scalar.h"
#include "../Ravioli/comp.h"
#include "../Ravioli/column.h"
#include "../LookUpTable/lookuptable.h"
#include "set.h"

////////////////
//            //
//  Property  //
//            //
////////////////
class Property : public Scalar {

  public:
    Property()                : con(true), cval(1.0) {} 
    Property(const char * nm) : con(true), cval(1.0) {nam=nm;} 
// ne treba: Property(const real & v) : con(true), mix(false), val(v)   {}

    real value() const;
    /* this is quite important: when declared as "virtual", it prevents
       calls to parent's "value" function, when "Property" is sent as a
       parameter to other functions (for plotting, for example) */
    virtual real value(const int i, const int j, const int k) const;
    void value(const real & v);
    void varies(const Domain & d);
    real value_comp(const int comp) const;

  protected:
    const Property * a;
    const Property * b;

  private:
    real   cval;
    bool   con;  
};

///////////////////
//               //
//  PropertyMix  // -> combination of two properties; interpolation
//               //
///////////////////
class PropertyMix : public Property {

  public:
    PropertyMix(const Property * pa, 
                const Property * pb, 
                const Scalar * coa, 
                const Scalar * cda,         
                const Scalar * cdb) {
      assert(pa !=NULL);
      assert(pb !=NULL);
      assert(coa!=NULL);
      a        = pa;
      b        = pb;
      c_a      = coa;
      c_disp_a = cda;
      c_disp_b = cdb;
    }

    real value(const int i, const int j, const int k) const;

  private:
    const Scalar  * c_a;
    const Scalar  * c_disp_a;
    const Scalar  * c_disp_b;
};

///////////////////
//               //
//  PropertyDiv  // -> combination of two properties; division
//               //
///////////////////
class PropertyDiv : public Property {

  public:
    PropertyDiv(const Property * pa, const Property * pb) {
      assert(pa !=NULL);
      assert(pb !=NULL);
      a  = pa;
      b  = pb;
    }

    real value(const int i, const int j, const int k) const {
      return a -> value(i,j,k) / b -> value(i,j,k);
    }
};

///////////////////
//               //
//  PropertyMul  // -> combination of two properties; multiplication
//               //
///////////////////
class PropertyMul : public Property {

  public:
    PropertyMul(const Property * pa, const Property * pb) {
      assert(pa !=NULL);
      assert(pb !=NULL);
      a  = pa;
      b  = pb;
    }

    real value(const int i, const int j, const int k) const {
      return a -> value(i,j,k) * b -> value(i,j,k);
    }
};

#endif
