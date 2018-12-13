#ifndef PHASECHANGEVOF_H
#define PHASECHANGEVOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"
//#include "nucleation.h"

//#define IB

class PhaseChangeVOF : public Centered {
  public:
    PhaseChangeVOF(const Scalar & mdot, 
                   const Scalar & mflx,
                   const Scalar & tpr,
                   const Scalar & tprs,
                   const Scalar & clr,
                   const Scalar & clrs,
                   const Scalar & vs,
                   const Vector & u, 
                   const Scalar & nx,
                   const Scalar & ny,
                   const Scalar & nz,
                   const Scalar & adens,
                   const Vector & fs,
                   const TIF & tifmodel,
                   Times & t,
                   Matter * flu,
                   real lat,
                   Matter * sol = NULL);
 
    ~PhaseChangeVOF();
    void update(const Scalar * diff_eddy = NULL);

    real get_turbP(){return turbP;}
    void set_turbP(real a){
      turbP=a;
      boil::oout<<"EnthalpyFD:turbP= "<<turbP<<"\n";
    }

    void cal_massflux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();
    void modify_vel(Vector & uvw, 
                    const Vector & cnew, const Vector & cold);

    bool get_upwind_flag() { return upwind_flag; }
    void set_upwind_flag(const bool flag = true) { upwind_flag = flag; }

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
    real mdot_cut(real m, real c);
    void mdot_cut();
    void mdot();
    void cal_norm_vect();
    void cal_gradt(const Scalar * diff_eddy = NULL);
    void gradt(const Scalar * diff_eddy = NULL);
    void gradt_ib(const Scalar * diff_eddy = NULL);

    void ext_gradt(Scalar & sca, const int iext, const Comp mcomp);
    void insert_bc_ext(const Comp mcomp);

    void setflag();

    void insert_bc_gradt(const Scalar * diff_eddy = NULL);
    void prepare_gradt8();

    real grad_2nd(const real tm0, const real tm1, const real tm2,
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

    void sources_clrs();
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

    bool upwind_flag;

    ScalarInt iflag;
    Scalar txv, tyv, tzv;
    Scalar txl, tyl, tzl;
    Scalar tpr, tprs, clr, clrs;
    Scalar adens;
    Scalar mx, my, mz;
    Scalar nx, ny, nz;
    Scalar stmp,stmp2,delta;
    Vector fs;
    Vector bndtpr;
    Scalar M;

    const TIF & tifmodel;

    real latent, rhol, rhov, lambdal, lambdav, cpl, cpv;
    real rhoave, rhodlt;
    real dxmin, clrsurf, pi;
    real turbP;
    real epsl; 
    real smdot_pos, smdot_neg, smdot_pos_macro, smdot_neg_macro;
};	

#endif

