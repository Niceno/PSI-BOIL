#ifndef PHASECHANGEVOF_H
#define PHASECHANGEVOF_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"
//#include "nucleation.h"

//#define IB

class PhaseChangeVOF : public Centered {
  public:
    PhaseChangeVOF(const Scalar & mdot, 
                   const Scalar & tpr,
                   const Scalar & tprs,
                   const Scalar & clr,
                   const Scalar & clrs,
                   const Scalar & vs,
                   const Vector & u, 
                   const Scalar & nx,
                   const Scalar & ny,
                   const Scalar & nz,
                   const Vector & fs,
                   Times & t,
                   Matter * flu,
                   real lat, real ts,
                   Scalar * tif = NULL,
                   Matter * sol = NULL);
 
    ~PhaseChangeVOF();
    void update(const Scalar * diff_eddy = NULL);

    real get_turbP(){return turbP;}
    void set_turbP(real a){
      turbP=a;
      boil::oout<<"EnthalpyFD:turbP= "<<turbP<<"\n";
    }

    Scalar M; // moved to public
    void cal_massflux(const Scalar * diff_eddy = NULL);
    void initialize();
    void finalize();
    void set_tif(Scalar * newtif) { tif = newtif; }
    void modify_vel(Vector & uvw, 
                    const Scalar & cnew, const Scalar & cold);

    void set_gradclr_ext(const Scalar & heavi);

  private:
    bool Interface(const int i, const int j, const int k,
                   const Comp m, const int dir);
#if 0
    bool Interface_x(const int i, const int j, const int k, const int dir);
    bool Interface_y(const int i, const int j, const int k, const int dir);
    bool Interface_z(const int i, const int j, const int k, const int dir);
    bool Interface1D_x(const int i, const int j, const int k, const int dir);
    bool Interface1D_y(const int i, const int j, const int k, const int dir);
    bool Interface1D_z(const int i, const int j, const int k, const int dir);
    bool Interface(const int i, const int j, const int k);
#endif
    void m(const Scalar * diff_eddy = NULL);
    real mdot_cut(real m, real c);
    void mdot_cut();
    void mdot();
    void cal_gradclr();
    void cal_norm_vect();
    void cal_gradt(const Scalar * diff_eddy = NULL);
    void ext_gradt(Scalar & sca, const int iext);
    void gradt (const Scalar * diff_eddy = NULL);

    void setflag();

    void prepare_gradt8();
    real grad_2nd(const real tm0, const real tm1, const real tm2,
                  const real dm1, const real dm2);
    real gradtx8(const int dir, const int i, const int j, const int k);
    real gradty8(const int dir, const int i, const int j, const int k);
    real gradtz8(const int dir, const int i, const int j, const int k);

    void sources_clrs();
    void sources_fext();
    void sources_sum();

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

    ScalarInt iflag;
    Scalar txv, tyv, tzv;
    Scalar txl, tyl, tzl;
    Scalar tpr, tprs, clr, clrs;
    Scalar gradclr;
    Scalar mx, my, mz;
    Scalar nx, ny, nz;
    Scalar stmp,stmp2,delta;
    Scalar * tif;
    Vector fs;
    Vector bndtpr;

    real tempnull; /* criterion for accepting bndtpr value as real */
    real latent, tsat, rhol, rhov, lambdal, lambdav, cpl, cpv;
    real rhoave, rhodlt;
    real dxmin, phisurf, pi;
    real turbP;
    real epsl; 
    real smdot_pos, smdot_neg, smdot_pos_macro, smdot_neg_macro;
};	

#endif

