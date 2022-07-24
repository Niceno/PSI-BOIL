#ifndef MICROLAYER_H
#define MICROLAYER_H

#include "../nucleation.h"
#include "../header.h"

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
                CommonHeatTransfer * cht,
                Heaviside * heavi,
                const Times * t,
                const real rs, 
                const real dmin, const real dmax = boil::unreal,
                Scalar * qsrc = NULL,
                const Sign sig = Sign::pos() );
    ~Microlayer() {};

    virtual void init() { 
      dmicro = boil::unreal;
      upkeep_after_seeding();
    }

    void update(real & smdot_micro,
                real & smdot_pos_macro_overwrite,
                real & smdot_neg_macro_overwrite);
    virtual void upkeep_after_seeding() { 
      upkeep_after_advance();
    }

    void upkeep_after_advance() {
#ifdef USE_VOF
      upkeep_after_advance_vof();
#else
      upkeep_after_advance_cip();
#endif
    }

    void hide_vf();

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

    /* heat transfer resistance */
    inline void set_hresis(const real h) {
      hresis = h;
      boil::oout<<"Microlayer::hresis= "<<hresis<<boil::endl;
      return;
    }
    inline real get_hresis() const {
      return hresis;
    }

    /* heat flux */
    inline real get_hflux_micro(const Dir d) const {return hflux_micro[int(d)];}
    inline real get_hflux_area(const Dir d) const {return area_sum[int(d)];}
    inline real get_hflux_area_l(const Dir d) const {return area_l[int(d)];}
    inline real get_hflux_area_v(const Dir d) const {return area_v[int(d)];}

    /* initial thickness */
    real d0(const int i, const int j, const int k) const;
    real d0max(const Comp m, const int i, const int j, const int k) const;

    /* clr & fs */
    void update_at_walls(Scalar & clr, Vector & fs);

  protected:
#ifdef USE_VOF
    void upkeep_after_advance_vof();
#else
    void upkeep_after_advance_cip();
    void store_dSprev();
    Scalar dSprev;
    bool str_dSprev;
#endif

    Scalar dmicro;
    Scalar * mdot, * tprs;

    real hresis;
    real slope, exp_slope;
    real rmax;
    real dmicro_min, dmicro_max;

    real * hflux_micro;
    real * area_sum, * area_l, * area_v;
};

#endif
