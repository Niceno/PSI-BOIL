#ifndef MICROLAYER_H
#define MICROLAYER_H

#include "../nucleation.h"
#include "../../Tifmodel/Schrage/schrage.h"

////////////////////////
//                    //
//  Microlayer model  //
//                    //
////////////////////////
class Microlayer : public Nucleation {
  public:
    Microlayer( Scalar & dmicro,
                Scalar * mdot,
                Scalar * tprs,
                const Scalar * tpr,
                Topology * topo,
                Heaviside * heavi,
                const TIF & tifmodel,
                const Times * t,
                Matter * f, const real rs, 
                const real dmin, const real dmax = boil::unreal,
                Matter * s = NULL,
                Scalar * qsrc = NULL,
                const Sign sig = Sign::pos() );
    ~Microlayer() {}

    virtual void init() { dmicro = boil::unreal; }

    void update(real & smdot_micro,
                real & smdot_pos_macro_overwrite,
                real & smdot_neg_macro_overwrite);
    virtual void upkeep_after_seeding();

    inline void set_slope(real r){ slope=r;
           boil::oout<<"Microlayer:slope is modified to "<<r<<"\n";};
    inline void set_slope(real r1, real r2) { slope=r1;
           exp_slope=r2;
           boil::oout<<"Microlayer:slope is modified to "<<r1<<"*r^"<<r2<<"\n";};
    inline real get_slope() const { return (slope); };
    inline real get_exp_slope() const { return (exp_slope); };

    inline void set_rmax(real r){ rmax = r;
            boil::oout<<"Microlayer:rmax is modified to "<<r<<"\n";};
    inline real get_rmax() const { return rmax;};

    /* heat flux */
    inline real get_hflux_micro(const Dir d) const {return hflux_micro[int(d)];}
    inline real get_hflux_area(const Dir d) const {return area_sum[int(d)];}
    inline real get_hflux_area_l(const Dir d) const {return area_l[int(d)];}
    inline real get_hflux_area_v(const Dir d) const {return area_v[int(d)];}

    inline const Matter * solid() const {return sld;}

  protected:
    void area_effect();
    void store_dSprev();
    real d0(const int i, const int j, const int k);

    Scalar dmicro;
    Scalar dSprev;
    Scalar * mdot, * tprs;
    const Domain * dom;
    const TIF * tifmodel;
    Matter * sld;

    real hresis;
    real slope, exp_slope;
    real rmax;
    real dmicro_min, dmicro_max;
    
    bool str_dSprev;

    real * hflux_micro;
    real * area_sum, * area_l, * area_v;
};

#endif
