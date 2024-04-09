#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Global/global_constants.h"
#include "../../SimulationTime/simulation_time.h"
#include "../../Matter/matter.h"
#include "site.h"

/////////////
//         //
//  Model  //
//         //
/////////////
class Nucleation {

  public:
    Nucleation ( Scalar * c, Scalar * tpr, Scalar * qsrc,
                 const Times * t, Scalar & dmicro,
                 Matter * f, const real rs, const real dm, 
                 const real l, const real ca );
    ~Nucleation();

    void save(const char *, const int);
    void load(const char *, const int);
    void rm  (const char *, const int);
    void save(std::ofstream &);
    void load(std::ifstream &);

    int  size() const { return  sites.size(); }
    int dsize() const { return dsites.size(); }

    void plant();
    void replant();
    void cutneck(const real r);
    void add(const Site & s);
    void st_active();
    real area_vapor(const int i, const int j, const int k, const Dir d);

    real clr_site  (const int i);
    real tpr_site  (const int i);
    real dmicro0(const int i, const int j, const int k );

    void set_optdat(bool b){ boptdat=b;
           boil::oout<<"nucleation:set_optdat = "<<b<<"\n";};
    real get_optdat(){ return (boptdat); };
    void optdat();

    void set_seed_period(real r){ seed_period=r; };
    real get_seed_period(){ return (seed_period); };

    void set_slope(real r){ slope=r; 
           boil::oout<<"nucleation:slope is modfied to "<<r<<"\n";};
    void set_slope(real r1, real r2) { slope=r1;
           slope=r2;
           b_slope=r2;
           boil::oout<<"nucleation:slope is modfied to "<<r1<<"*r+"<<r2<<"\n";};
    void set_exp_slope(real r1, real r2) { slope=r1;
           exp_slope=r2;
           boil::oout<<"nucleation:slope is modfied to "<<r1<<"*r^"<<r2<<"\n";};
    real get_slope(){ return (slope); };
    real get_b_slope(){ return (b_slope); };
    real get_exp_slope(){ return (exp_slope); };

    void set_rmax(real r){ rmax = r;
            boil::oout<<"nucleation:rmax is modified to "<<r<<"\n";};
    real get_rmax(){ return rmax;};

    void set_range_zoning(real r){range_zoning = r;
            boil::oout<<"nucleation:range_zoning is modified to "<<r<<"\n";};
    real get_range_zoning(){return (range_zoning);};

    void set_heat_sink(bool b){heat_sink = b;
            boil::oout<<"nucleation:heat_sink is modified to "<<b<<"\n";};
    real get_heat_sink(){return (heat_sink);};

    real cangle() {return cang;};
    real sigma()  {return sgm;}

    std::vector<Site> sites, dsites;
    Scalar dmicro;
    real ** dSprev;
    bool store_dSprev;
    real dmicro_min;

    const Matter * fluid() const {return flu;}

  protected:
    void set_range(std::vector<Site> & s);
    real zftmin    (Range<real> xr, Range<real> yr, Range<real> zr);
    real area_vapor_sum(Range<real> xr, Range<real> yr, Range<real> zr);
    void zoning();

    Scalar * clr;
    Scalar * tpr;
    Scalar * qsrc;
    Matter * flu;
    const Times * time;
    real seed_period; //active_tpr, zplant;
    real period_cut_replant;  // period of preventing replant after cutneck
    real dxmin, zbtm;
    real cang, sgm;  // contact angle, surface tension coefficient

  private:
    void dummy_add(const Site & s){ dsites.push_back(s); 
                              set_range(dsites);};
    real slope, b_slope, exp_slope, latent, rseed;
    real rmax;
    bool bzoning,boptdat,heat_sink;
    real clrsurf, eps_clr;
    std::vector<int> id_nearRegion, idd_nearRegion;
    std::ofstream ofstpr,ofsclr;
    real range_zoning;
};

#endif
