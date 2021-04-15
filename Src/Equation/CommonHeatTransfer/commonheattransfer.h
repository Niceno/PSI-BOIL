#ifndef CHTT_H
#define CHTT_H

#include "../../Global/global_func.h"
#include "../Topology/topology.h"
#include "../Tifmodel/tif.h"

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
                       Matter * flu, Matter * sol = NULL);

    ~CommonHeatTransfer();

#include "commonheattransfer_inline.h"

    /* cell distances */
    real distance_center(const Sign sig, const Comp & m,
                         const int i, const int j, const int k) const;
    real distance_face(const Sign sig, const Comp & m,
                       const int i, const int j, const int k) const;

    /* distances to interface */
    real distance_int(const Sign dir, const Comp & m,
                      const int i, const int j, const int k,
                      real & tint, Sign & cell_marker,
                      const ResistEval re, const Old old) const;
    real distance_int_x(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint, Sign & cell_marker,
                        const ResistEval re, const Old old) const;
    real distance_int_y(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint, Sign & cell_marker,
                        const ResistEval re, const Old old) const;
    real distance_int_z(const Sign dir,
                        const int i, const int j, const int k,
                        real & tint, Sign & cell_marker,
                        const ResistEval re, const Old old) const;

    /* inverse interfacial heat transfer resistance */
    real evaluate_resinv(const Sign dir, const Comp & m,
                         const int i0, const int j0, const int k0,
                         const int ii, const int ji, const int ki,
                         const real dist) const;

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

    /* 
    real gradt_ib(const Sign dir, const Comp & mcomp,
                  const int i, const int j, const int k,
                  const Old old,
                  Scalar & val) const; */

    /* calculating a variable-stencil difference */
    real first_derivative(const bool is_solid, const Comp & m,
                          const int i, const int j, const int k,
                          const AccuracyOrder & accuracy_order,
                          const bool discard_points,
                          const Old old = Old::no) const;
    real second_derivative(const bool is_solid, const Comp & m,
                           const int i, const int j, const int k,
                           const AccuracyOrder & accuracy_order,
                           const bool discard_points,
                           const Old old = Old::no) const;

    void construct_stencil(std::vector<StencilPoint> & stencil,
                           const bool is_solid, const Comp & m,
                           const int i, const int j, const int k,
                           const AccuracyOrder & accuracy_order,
                           const bool discard_points,
                           const Old old) const;

    /* calculate solid wall temperature */
    void calculate_node_temperature(const Scalar * diff_eddy = NULL);

    Topology * topo;
    const TIF & tifmodel;

    const Matter * fluid() const {return flu;}
    const Matter * solid() const {return sol;}

    const Scalar & tmp() const {return tpr;}
    const Vector & node_tmp_sol() const {return bndtpr_sol;}
    const Vector & node_tmp_flu() const {return bndtpr_flu;}
    Vector & node_tmp_sol() {return bndtpr_sol;}
    Vector & node_tmp_flu() {return bndtpr_flu;}

  private:
    /* effect of interfacial heat transfer resistance */
    void resTint(const Sign & dir, const Comp & m,
                 const int i0, const int j0, const int k0,
                 const int ii, const int ji, const int ki,
                 const int i1, const int j1, const int k1,
                 const real dist, const Sign & cell_marker,
                 real & tint, const Old old) const;

    /* one point for difference stencil construction */
    bool add_point(const int i0, const int j0, const int k0,
                   const int i1, const int j1, const int k1,
                   const Sign dir, const Comp & m,
                   const bool is_solid, bool & terminate,
                   bool & interface_reached,
                   std::vector<StencilPoint> & stencil,
                   const Old old) const;


    Scalar tpr;
    Vector bndtpr_sol, bndtpr_flu;

    Matter * flu;
    Matter * sol;

    boil::func_ijk_real dirac_wall_source_func;

    real val_rhov,val_rhol,val_cpv,val_cpl,val_lambdav,val_lambdal;
    real turbP; /* turbulent Prandtl number */
    real dirac_wall_source_val; /* units W/m2 */ 
    real wall_resistance_val; /* solid-fluid contact resistance */
    /* heat transfer resistance, only for liquid */
    real int_resistance_liq_val; 
    bool use_int_resist;
    bool dirac_wall_source_variable;
};

#endif
