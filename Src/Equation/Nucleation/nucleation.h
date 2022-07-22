#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Global/global_constants.h"
#include "../../Global/global_realistic.h"
#include "../../SimulationTime/simulation_time.h"
#include "../../Matter/matter.h"
#include "../CommonHeatTransfer/commonheattransfer.h"
#include "../Heaviside/heaviside.h"
#include "site.h"
#include "header.h"

/////////////
//         //
//  Model  //
//         //
/////////////
class Nucleation {

  public:
    Nucleation ( CommonHeatTransfer * cht,
                 Heaviside * heavi,
                 const Times * t,
                 const real rs,
                 Scalar * qsrc = NULL, 
                 const Sign sig = Sign::pos() );
    /* sign positive -> liquid = 1 */
    ~Nucleation();

    void save(const char *, const int);
    void load(const char *, const int);
    void rm  (const char *, const int);
    void save(std::ofstream &);
    void load(std::ifstream &);

    inline int  size() const { return  sites.size(); }
    inline int dsize() const { return dsites.size(); }

    void plant();
    void replant();
    virtual void init() {}

#ifndef USE_VOF
    void cutneck(const real r);
    inline void set_cutneck_mult(real cm){ rcut = rseed*cm; };
    inline real get_cutneck_mult() const { return rcut; };
#endif

    void add(const Site & s);
    void st_active();
    virtual void upkeep_after_seeding() {};

    real clr_site  (const int i);
    real tpr_site  (const int i);

    real distance_from_site(const int i, const int j, const int k) const;

    //real area_vapor(const Dir d,
    //                const int i, const int j, const int k) const;
    real area_vapor(const Sign sig, const Comp & mcomp,
                    const int i, const int j, const int k) const;

    inline void set_seed_period(real r){ seed_period=r; };
    inline real get_seed_period() const { return (seed_period); };

    inline void set_prevent_replant_period(real r){ period_prevent_replant=r; };
    inline real get_prevent_replant_period() const { return (period_prevent_replant); };

    inline void set_threshold_c(real c){ threshold_c = c; };
    inline real get_threshold_c() const { return threshold_c; };

    inline void set_zoning_limiting(bool zl, real zlmult = 0.0) {
      limit_zoning = zl;
      zoning_limit_multiplier = zlmult;
    }
    inline bool get_zoning_limiting(real & zlmult) const {
      zlmult = zoning_limit_multiplier;
      return limit_zoning;
    }

    bool in_vapor(const int i, const int j, const int k) const;
    bool in_vapor(const real c) const;
    bool below_threshold(const int i, const int j, const int k) const;
    bool below_threshold(const real c) const;

    std::vector<Site> sites, dsites;

    //const Matter * fluid() const {return flu;}
    //const Topology * topology() const {return topo;}

  protected:
    void set_range(std::vector<Site> & s);

    real area_vapor_sum(Range<real> xr, Range<real> yr, Range<real> zr);

    void zoning();

    void dummy_add(const Site & s){ dsites.push_back(s); 
                              set_range(dsites);};
    real stratified_sphere(const int i, const int j, const int k,
                        const real xcent, const real ycent, const real zcent);

    void plant_site(const int ns, const bool seed_source = true);
    void plant_dummy_site(const int nsd);

    //Topology * topo;
    CommonHeatTransfer * cht;
    Heaviside * heavi;

    Scalar * vf;
    Scalar * clr;
    //const Scalar * tpr;
    Scalar * qsrc;
    //Matter * flu;
    const Times * time;

    real seed_period;
    real period_prevent_replant;  // period of preventing replant after cutneck
    real dxmin, eps;
    real zbtm;

    real rhol, rhov, lambdal, lambdav, cpl, cpv, latent, mmass;

    real rseed;
    bool bzoning;
    std::vector<int> id_nearRegion, idd_nearRegion;
    Sign matter_sig;

    bool limit_zoning;
    real zoning_limit_multiplier;
    real threshold_c;

#ifndef USE_VOF
    real rcut;
#endif
};

#endif
