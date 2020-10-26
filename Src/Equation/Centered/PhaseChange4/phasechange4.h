#ifndef PHASECHANGE4_H
#define PHASECHANGE4_H

#include <cmath>
#include <set>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"
#include "../../Topology/topology.h"
#include "../../../Ravioli/htwallmodel.h"

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
                 Topology * topo,
                 const TIF & tifmodel,
                 Times & t,
                 Matter * flu,
                 Matter * sol = NULL,
                 HTWallModel * htwallmodel = NULL,
                 Sign matter_sig = Sign::pos());
 
    ~PhaseChange4();
    void update(const Scalar * diff_eddy = NULL);

    void mass_flux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();
    inline void sources() { sources_vfs(); sources_fext(); sources_sum(); }

    real lambda(const int i, const int j, const int k,
                const Scalar * diff_eddy = NULL) const;
    real lambda_inv(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const;

#include "phasechange4_inline.h"

    Vector & node_tmp() {return bndtpr;}

    void request_flux(const int i, const int j, const int k,
                      std::vector<real> & tv,
                      std::vector<real> & tl) const {
      tv = {txv[i][j][k],tyv[i][j][k],tzv[i][j][k],tnv[i][j][k]};
      tl = {txl[i][j][k],tyl[i][j][k],tzl[i][j][k],tnl[i][j][k]};
      return;
    }


  private:
    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool interface(const int i, const int j, const int k);
    void m();
    real mdot_cut(real m, real c, real & mcut);
    real vfs_cut(real vfsval, real vfval);
    void mdot();

    void heat_flux(const Scalar * diff_eddy = NULL);
    void cal_hf(const Scalar * diff_eddy = NULL);
    void insert_bc_hf(const Scalar * diff_eddy);
    void dirac_source_terms();
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

    /* ghost distance */
    real ghost_distance(const Comp & m, const Sign & cell_marker,
                        const int i, const int j, const int k);
 
    real Tint(const int i, const int j, const int k);
 
    real gradt1D(const bool is_solid, const Comp & m,
                 const int i, const int j, const int k);

    bool add_point(const int i0, const int j0, const int k0,
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

    ScalarInt tempflag,iflag;
    Scalar txv, tyv, tzv;
    Scalar txl, tyl, tzl;
    Scalar tnv, tnl;

    Scalar tpr, tprs, vf, vfs;
    Scalar adens;
    Scalar nx;
    Scalar ny;
    Scalar nz;
    Vector bndtpr;
    Scalar M;

    const TIF & tifmodel;
    const Sign matter_sig; /* pos: liquid is phi=1 and vice versa */
    Topology * topo;
    HTWallModel * htwallmodel;
 
    real rhol, rhov, lambdal, lambdav, cpl, cpv;
    real turbP;
    real smdot_pos, smdot_neg;
    
    int accuracy_order;
    bool default_value_for_htwallmodel;
    bool use_unconditional_extrapolation,
         discard_points_near_interface;
};	

#endif

