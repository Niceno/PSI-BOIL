#ifndef MATTER_H
#define MATTER_H

#include <vector>

#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Ravioli/comp.h"
#include "../Ravioli/column.h"
#include "../LookUpTable/lookuptable.h"
#include "set.h"
#include "property.h"

//////////////
//          //
//  Matter  //
//          //
//////////////
class Matter {

  public:
/*  Matter() do not use this one, dangerous */        

    Matter(const Domain & d, const char * nm = NULL);

    //Matter(const Domain & d, const char * nm);

    Matter(const Matter & a, 
           const Matter & b, 
           const Scalar * ca,          /* concentration of component "a" */
           const Scalar * cda = NULL,  /* concentration of dispersed "a" */
           const Scalar * cdb = NULL); /* concentration of dispersed "b" */

    Matter(const Matter & a, 
           const Matter & b, 
           const Scalar * ca,          /* concentration of component "a" */
           const Vector * bdca,        /* concentration of "a", staggered */
           const Scalar * cda = NULL,  /* concentration of dispersed "a" */
           const Scalar * cdb = NULL); /* concentration of dispersed "b" */

    real rho   (const int i,
                const int j,
                const int k) const {return dens->value(i,j,k);}
    real rho   (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return dens->value(m,i,j,k);
                }

    real rho   (const int comp) const {return dens->value_comp(comp);}

    real mu    (const int i,
                const int j,
                const int k) const {return visc->value(i,j,k);}
    real mu    (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return visc->value(m,i,j,k);
                }
    real mu    (const int comp) const {return visc->value_comp(comp);}

    real cp    (const int i,
                const int j,
                const int k) const {return capa->value(i,j,k);}
    real cp    (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return capa->value(m,i,j,k);
                }
    real cp    (const int comp) const {return capa->value_comp(comp);}

    real lambda(const int i,
                const int j,
                const int k) const {return cond->value(i,j,k);}
    real lambda(const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return cond->value(m,i,j,k);
                }
    real lambda(const int comp) const {return cond->value_comp(comp);}

    real gamma (const int i,
                const int j,
                const int k) const {return diff->value(i,j,k);}
    real gamma (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return diff->value(m,i,j,k);
                }
    real gamma (const int comp) const {return diff->value_comp(comp);}

    real beta  (const int i,
                const int j,
                const int k) const {return texp->value(i,j,k);}
    real beta  (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return texp->value(m,i,j,k);
                }
    real beta  (const int comp) const {return texp->value_comp(comp);}

    real mmass (const int i,
                const int j,
                const int k) const {return molm->value(i,j,k);}
    real mmass (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                    return molm->value(m,i,j,k);
                }
    real mmass (const int comp) const {return molm->value_comp(comp);}


    real sigma (const int i,
                const int j,
                const int k) const {assert(tens); return tens->value(i,j,k);}
    real sigma (const Comp & m,
                const int i,
                const int j,
                const int k) const {
                  assert(tens);
                  return tens->value(m,i,j,k);
                }
    
    real latent(const int i,
                const int j,
                const int k) const {assert(heat); return heat->value(i,j,k);}
    real latent(const Comp & m,
                const int i,
                const int j,
                const int k) const {
                  assert(heat);
                  return heat->value(m,i,j,k);
                }

    /* set (initialize) physical properties */
    void rho   (const real & v) {dens->value(v);}
    void mu    (const real & v) {visc->value(v);}
    void cp    (const real & v) {capa->value(v);}
    void lambda(const real & v) {cond->value(v);}
    void gamma (const real & v) {diff->value(v);}
    void beta  (const real & v) {texp->value(v);}
    void mmass (const real & v) {molm->value(v);}

    void sigma (const real & v) {
      if(tens == NULL) {
        boil::oout << "# Fatal: specifying surface tension ";
        boil::oout << "makes sense only for mixtures. Exiting!"; 
        boil::oout << boil::endl;    
        exit(0);
      }
      tens->value(v);
    }
    void latent(const real & v) {
      if(heat == NULL) {
        boil::oout << "# Fatal: specifying latent heat ";
        boil::oout << "makes sense only for mixtures. Exiting!";
        boil::oout << boil::endl;
        exit(0);
      }
      heat->value(v);
    }

    
    const Property * rho()    const {return dens;}
    const Property * mu()     const {return visc;}
    const Property * cp()     const {return capa;}
    const Property * lambda() const {return cond;}
    const Property * gamma()  const {return diff;}
    const Property * beta()   const {return texp;}
    const Property * mmass()  const {return molm;}

    const Property * sigma()  const {return tens;}
    const Property * latent() const {return heat;}

    /* rescale properties by factors of length and time */
    void rescale(const real xmult = 1., const real tmult = 1.,
                 const real mmult = 1.);

    /* set certain property to be variable */
    void variable(const Set & s);

    /* set all properties to be variable */
    void variable();

    /* look up a property */
    void look_up(const Set & s, const Scalar & sca,
                 const LookUpTable & pt,  
                 const Column & col0, const Column & col1);

    /* check if matter is mixture */
    inline bool mixture() const { return mixt; }

//  friend std::ostream & 
//    operator << (std::ostream & os, const Property & prop);

  private:
    Property * dens;      /* [kg/m^3] */
    Property * visc;
    Property * capa;
    Property * cond;
    Property * diff;
    Property * texp;
    Property * molm;
    Property * tens;
    Property * heat;

    Property * dens_o_visc;

    const Domain * dom;

    std::string nam;

    bool mixt;
};

#endif
