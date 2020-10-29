#ifndef CHTT_H
#define CHTT_H

#include "../Topology/topology.h"
#include "../Tifmodel/tif.h"
#include "../../Ravioli/htwallmodel.h"

///////////////////////////
//                       //
// Common Heat Transfer  //
//                       //
///////////////////////////
/* 
 * this class encapsulates common methods of phase change and enthalpy
 */
class CommonHeatTransfer {
  public:
    CommonHeatTransfer(const Scalar & tpr,
                       Topology * topo, const TIF & tifmodel, 
                       Matter * flu, Matter * sol = NULL,
                       HTWallModel * htwallmodel = NULL);

    ~CommonHeatTransfer();

#include "commonheattransfer_inline.h"

    /* cell distances */
    real distance_center(const Sign sig, const Comp & m,
                         const int i, const int j, const int k) const;
    real distance_face(const Sign sig, const Comp & m,
                       const int i, const int j, const int k) const;

    /* new distances to interface */
    inline real distance_int(const Sign dir, const Comp & m,
                             const int i, const int j, const int k,
                             real & tint, const Old old) const {
      return (old==Old::yes) ?
             distance_int_old(dir,m,i,j,k,tint) :
             distance_int(dir,m,i,j,k,tint);
    };
    real distance_int(const Sign dir, const Comp & m,
                      const int i, const int j, const int k,
                      real & tint) const;
    real distance_int_x(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint) const;
    real distance_int_y(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint) const;
    real distance_int_z(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint) const;

    /* old distances to interface */
    real distance_int_old(const Sign dir, const Comp & m,
                          const int i, const int j, const int k,
                          real & tint) const;
    real distance_int_x_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint) const;
    real distance_int_y_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint) const;
    real distance_int_z_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint) const;

    /* ghost distance */
    real ghost_distance(const Comp & m, const Sign & cell_marker,
                        const int i, const int j, const int k) const;

    /* test domain edge */
    bool edge(const Sign dir, const Comp & m,
              const int i, const int j, const int k) const;

    /* thermal conductivity */
    real lambda(const int i, const int j, const int k,
                const Scalar * diff_eddy = NULL) const;
    real lambda_inv(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const;
    real lambda_old(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const;
    real lambda_inv_old(const int i, const int j, const int k,
                        const Scalar * diff_eddy = NULL) const;

    /* hflux and gradient */
    real hflux_wall(const Scalar & s, const Dir d,
                    const Scalar * diff_eddy = NULL) const;
    real hflux_wall_ib(const Scalar & s, 
                       const Scalar * diff_eddy = NULL) const;

    real gradt_ib(const Sign dir, const Comp & mcomp,
                  const int i, const int j, const int k,
                  const Old old,
                  Scalar & val) const;

    /* calculating a variable-stencil gradient */
    real gradt1D(const bool is_solid, const Comp & m,
                 const int i, const int j, const int k,
                 const AccuracyOrder & accuracy_order,
                 const bool discard_points,
                 const Old old = Old::no) const;
    bool add_point(const int i0, const int j0, const int k0,
                   const int i1, const int j1, const int k1,
                   const Sign dir, const Comp & m,
                   const bool is_solid, bool & terminate,
                   std::vector<real> & stencil,
                   std::vector<real> & values,
                   const Old old) const;

    /* calculate solid wall temperature */
    void calculate_node_temperature(const Scalar * diff_eddy = NULL);

    /* members */
    inline HTWallModel & heat_transfer_wall_model() {
      return *htwallmodel;
    }
    inline const HTWallModel & heat_transfer_wall_model() const {
      return *htwallmodel;
    }
    Topology * topo;
    const TIF & tifmodel;

    const Matter * fluid() const {return flu;}
    const Matter * solid() const {return sol;}

    const Scalar & tmp() const {return tpr;}
    const Vector & node_tmp() const {return bndtpr;}
    Vector & node_tmp() {return bndtpr;}

  private:
    Scalar tpr;
    Vector bndtpr;

    Matter * flu;
    Matter * sol;
    HTWallModel * htwallmodel;

    bool default_value_for_htwallmodel;

    real val_rhov,val_rhol,val_cpv,val_cpl,val_lambdav,val_lambdal;
    real turbP; /* turbulent Prandtl number */
    real resistance_equivalent; /* length scale of heat transfer resistance */
};

#endif
