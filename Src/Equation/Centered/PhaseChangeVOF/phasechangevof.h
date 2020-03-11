#ifndef PHASECHANGEVOF_H
#define PHASECHANGEVOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"
#include "../../Topology/topology.h"
//#include "nucleation.h"

#define IB

class PhaseChangeVOF : public Centered {
  public:
    PhaseChangeVOF(const Scalar & mdot, 
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
 
    ~PhaseChangeVOF();
    void update(const Scalar * diff_eddy = NULL);

    void cal_massflux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();
    void modify_vel(Vector & uvw, 
                    const Vector & cnew, const Vector & cold);
    void modify_vel(Vector & uvw, 
                    const Scalar & cnew, const Scalar & cold);

    inline real get_turbP() const { return turbP; }
    inline void set_turbP(const real a) {
      turbP = a;
      boil::oout<<"PhaseChangeVOF::turbP= "<<turbP<<"\n";
    }

    inline bool get_upwind_flag() const { return upwind_flag; }
    inline void set_upwind_flag(const bool flag = true) {
      upwind_flag = flag;
      boil::oout<<"PhaseChangeVOF::upwind_flag= "<<upwind_flag<<"\n";
    }

    inline bool get_near_wall_modelling() const { return near_wall_modelling; }
    inline void set_near_wall_modelling(const bool flag) {
      near_wall_modelling = flag; 
      boil::oout<<"PhaseChangeVOF::near_wall_modelling= "<<near_wall_modelling
                <<"\n";
    }

  private:
    bool Interface(const int dir, const Comp m,
                   const int i, const int j, const int k);
    bool Interface(const int i, const int j, const int k);
#if 0
    bool Interface_x(const int i, const int j, const int k, const int dir);
    bool Interface_y(const int i, const int j, const int k, const int dir);
    bool Interface_z(const int i, const int j, const int k, const int dir);
    bool Interface1D_x(const int i, const int j, const int k, const int dir);
    bool Interface1D_y(const int i, const int j, const int k, const int dir);
    bool Interface1D_z(const int i, const int j, const int k, const int dir);
#endif
    void m(const Scalar * diff_eddy = NULL);
    real mdot_cut(real m, real c, real & mcut);
    void mdot();

    void cal_gradt_fluid(const Scalar * diff_eddy = NULL);
    void cal_gradt_ib(const Scalar * diff_eddy = NULL);
    void correct_gradt_at_ib(const Scalar * diff_eddy = NULL);
    void insert_bc_gradt_at_walls(const Scalar * diff_eddy = NULL);
    void gradt(const Scalar * diff_eddy = NULL);

    void near_wall_model(const Scalar * diff_eddy = NULL);

    void calculate_node_temperature(const Scalar * diff_eddy = NULL);

    real second_order_difference(const real tm0, const real tm1, const real tm2,
                                 const real dm1, const real dm2);
    real gradtx(const int dir, const int i, const int j, const int k);
    real gradty(const int dir, const int i, const int j, const int k);
    real gradtz(const int dir, const int i, const int j, const int k);
    real gradtx8(const int dir, const int i, const int j, const int k);
    real gradty8(const int dir, const int i, const int j, const int k);
    real gradtz8(const int dir, const int i, const int j, const int k);
    real gradtx9(const int dir, const int i, const int j, const int k);
    real gradty9(const int dir, const int i, const int j, const int k);
    real gradtz9(const int dir, const int i, const int j, const int k);

    void sources_vfs();
    void sources_fext();
    void sources_sum();

    real distance(const int i, const int j, const int k,
                  const int dir, const Comp m, real & tint);
    real distance_x(const int i, const int j, const int k,
                    const int dir, real & tint);
    real distance_y(const int i, const int j, const int k,
                    const int dir, real & tint);
    real distance_z(const int i, const int j, const int k,
                    const int dir, real & tint);
    bool distance1D_x(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
    bool distance1D_y(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
    bool distance1D_z(const int i, const int j, const int k,
                      const int dir, real & tint, real & dist);
 
    real Tint(const int i, const int j, const int k);

    real temperature_node(real len_s, real lam_s, real tmp_s
                        , real len_f, real lam_f, real tmp_f);
 
    real lambda(const int i, const int j, const int k,
                const Scalar * diff_eddy = NULL) const;

    bool upwind_flag, near_wall_modelling;

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
    real smdot_pos, smdot_neg, smdot_pos_macro, smdot_neg_macro;
};	

#endif

