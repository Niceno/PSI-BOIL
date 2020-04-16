#ifndef PHASECHANGE4_H
#define PHASECHANGE4_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"
#include "../../Topology/topology.h"

#define IB

class PhaseChange4 : public Centered {
  public:
    PhaseChange4(const Scalar & mdot, 
                 const Scalar & mflx,
                 const Scalar & tpr,
                 const Scalar & tprs,
                 const Scalar & vf,
                 const Scalar & vfs,
                 const Scalar & vs,
                 const Vector & u, 
                 Topology & topo,
                 const TIF & tifmodel,
                 Times & t,
                 Matter * flu,
                 Matter * sol = NULL,
                 Sign sig = Sign::pos());
 
    ~PhaseChange4();
    void update(const Scalar * diff_eddy = NULL);

    void mass_flux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();

    inline real get_turbP() const { return turbP; }
    inline void set_turbP(const real a) {
      turbP = a;
      boil::oout<<"PhaseChange4::turbP= "<<turbP<<"\n";
    }

    /* testing */
    bool test_differences(const int count);
    bool test_differences(const std::vector<real> & stencil,
                          const std::vector<real> & coefficients);
    real evaluate_polynomial(const int order,
                             const std::vector<real> & coefficients,
                             const real x);


  private:
    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool interface(const int i, const int j, const int k);
    void m();
    real mdot_cut(real m, real c, real & mcut);
    void mdot();

    void heat_flux(const Scalar * diff_eddy = NULL);
    void cal_hf(const Scalar * diff_eddy = NULL);
    void calculate_node_temperature(const Scalar * diff_eddy = NULL);

    void sources_vfs();
    void sources_fext();
    void sources_sum();

    real distance_int(const Sign dir, const Comp & m,
                      const int i, const int j, const int k,
                      real & tint);
    real distance_int_x(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint);
    real distance_int_y(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint);
    real distance_int_z(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint);
    bool distance1D_int_x(const int i, const int j, const int k,
                          const Sign dir, real & tint, real & dist);
    bool distance1D_int_y(const int i, const int j, const int k,
                          const Sign dir, real & tint, real & dist);
    bool distance1D_int_z(const int i, const int j, const int k,
                          const Sign dir, real & tint, real & dist);
 
    real Tint(const int i, const int j, const int k);

    real temperature_node(real len_s, real lam_s, real tmp_s,
                          real len_f, real lam_f, real tmp_f);
 
    real lambda(const int i, const int j, const int k,
                const Scalar * diff_eddy = NULL) const;
    real lambda_inv(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const;

    real gradt1D(const bool is_solid, const Comp & m,
                 const int i, const int j, const int k);

    void add_point(const int i0, const int j0, const int k0,
                   const int i1, const int j1, const int k1,
                   const Sign dir, const Comp & m,
                   const bool is_solid, bool & terminate,
                   std::vector<real> & stencil,
                   std::vector<real> & values);

    real distance_center(const Sign sig, const Comp & m,
                         const int i, const int j, const int k);
    real distance_face(const Sign sig, const Comp & m,
                       const int i, const int j, const int k);
    bool edge(const Sign dir, const Comp & m,
              const int i, const int j, const int k);

    real second_order_difference(const std::vector<real> & stencil,
                                const std::vector<real> & values);
    real third_order_difference(const std::vector<real> & stencil,
                                const std::vector<real> & values);
    real fourth_order_difference(const std::vector<real> & stencil,
                                 const std::vector<real> & values);

    ScalarInt tempflag,iflag;
    Scalar txv, tyv, tzv;
    Scalar txl, tyl, tzl;
    Scalar tnv, tnl;

    Scalar tpr, tprs, vf, vfs;
    Scalar clr,adens;
    Scalar nx;
    Scalar ny;
    Scalar nz;
    Vector fs;
    Vector bndtpr;
    Scalar M;

    const TIF & tifmodel;
    const Sign sig; /* pos: liquid is phi=1 and vice versa */
    Topology * topo;
 
    real rhol, rhov, lambdal, lambdav, cpl, cpv;
    real clrsurf;
    real turbP;
    real epsl; 
    real smdot_pos, smdot_neg;
};	

#endif

