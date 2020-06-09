#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Global/global_constants.h"
#include "../../Global/global_realistic.h"
#include "../../SimulationTime/simulation_time.h"
#include "../../Matter/matter.h"
#include "../Topology/topology.h"
#include "../Heaviside/heaviside.h"
#include "site.h"

/////////////
//         //
//  Model  //
//         //
/////////////
class Nucleation {

  public:
    Nucleation ( Topology * topo,
                 Heaviside * heavi,
                 const Scalar * tpr,
                 const Times * t,
                 Matter * f, const real rs,
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
    void cutneck(const real r);
    void add(const Site & s);
    void st_active();

    real clr_site  (const int i);
    real tpr_site  (const int i);

    real distance_from_site(const int i, const int j, const int k) const;

    //real area_vapor(const Dir d,
    //                const int i, const int j, const int k) const;
    real area_vapor(const Sign sig, const Comp & mcomp,
                    const int i, const int j, const int k) const;

    inline void set_seed_period(real r){ seed_period=r; };
    inline real get_seed_period() const { return (seed_period); };

    inline void set_cutneck_mult(real cm){ rcut = rseed*cm; };
    inline real get_cutneck_mult() const { return rcut; };

    bool in_vapor(const int i, const int j, const int k) const;
    bool in_vapor(const real c) const;

    std::vector<Site> sites, dsites;

    const Matter * fluid() const {return flu;}
    const Topology * topology() const {return topo;}

  protected:
    void set_range(std::vector<Site> & s);

    real area_vapor_sum(Range<real> xr, Range<real> yr, Range<real> zr);

    void zoning();

    void dummy_add(const Site & s){ dsites.push_back(s); 
                              set_range(dsites);};
    real stratified_sphere(const int i, const int j, const int k,
                        const real xcent, const real ycent, const real zcent);

    virtual void upkeep_after_seeding() {};
    void plant_site(const int ns, const bool seed_source = true);
    void plant_dummy_site(const int nsd);

    Topology * topo;
    Heaviside * heavi;

    Scalar * vf;
    const Scalar * clr;
    const Scalar * tpr;
    Scalar * qsrc;
    Matter * flu;
    const Times * time;

    real seed_period;
    real period_cut_replant;  // period of preventing replant after cutneck
    real dxmin, eps;
    real zbtm;

    real rhol, rhov, lambdal, lambdav, latent, mmass;

    real rseed, rcut;
    bool bzoning;
    std::vector<int> id_nearRegion, idd_nearRegion;
    Sign sig;
};

#endif
