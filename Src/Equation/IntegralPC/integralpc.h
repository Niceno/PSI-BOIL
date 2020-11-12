#ifndef INTEGRALPC_H
#define INTEGRALPC_H

#include "../Staggered/Momentum/momentum.h"
#include "../Centered/EnthalpyFD/enthalpyfd.h"
#include "../Centered/PhaseChange4/phasechange4.h"
#include "../Centered/CavityPressure/cavitypressure.h"
#include "../../Solver/Additive/additive.h"

/////////////////////////////
//                         //
//  Integral Phase Change  //
//                         //
/////////////////////////////
/* A class for solving combined phase-change with inertial effects */
class IntegralPC {
  public:
    IntegralPC(Momentum & NS,
               EnthalpyFD & EFD, PhaseChange4 & PC, 
               CavityPressure & CAPR, AC & MC,
               TIF & TSAT, const Times & T,
               Scalar & F, Scalar & PRESS,
               Scalar & TPR, const Scalar & TPR_OLD,
               Vector & UVW, const Vector & UVW_OLD,
               Vector & UVW_CAV,
               const std::function<real(const real t)> & MODEL) :
      ns(NS), enthFD(EFD), pc(PC), capr(CAPR), multigrid_cavity(MC),
      tsat(TSAT), time(T),
      f(F), press(PRESS),
      tpr(TPR), tpr_old(TPR_OLD), 
      uvw(UVW), uvw_old(UVW_OLD), uvw_cav(UVW_CAV), model(MODEL)
    {
      c0 = Cycle::none();
      c1 = Cycle::F();
      rr_cav = ResRat(5e-5);
      rt_cav = ResTol(5e-5);
      
      const int niter = 1000;
      MaxIter mm = MaxIter(niter);
      mi = {mm,mm,mm};

      rr_ns = ResRat(1e-14);
      
      errtol = ResTol(1e-6);
    }

    ~IntegralPC() {};

    real solve(const ResTol & rt, const real pinf, const Range<real> & tprr,
               const bool progress = true);

    void set_mg_params(const Cycle & C0,  const Cycle & C1, const ResTol & RT, 
                       const ResRat & RR, const std::array<MaxIter,3> & MI) {
      c0 = C0;
      c1 = C1;
      rr_cav = RR;
      rt_cav = RT;
      mi = MI;
      
      return;
    }

    void set_ns_params(const ResRat & RR) {
      rr_ns = RR; return;
    }

    void set_errtol(const ResTol & RT) {
      errtol = RT; return;
    }

    void set_temperature(const real T) {
      tsat.set_tref(T);
      t_old = T; return;
    }

    real get_temperature() const {
      return t_old;
    }

  private:
    real cavity_evaluation(const real pinf, const real pcc);
    void momentum_evaluation();
    real phase_change_evaluation(const real tval);

    real evaluate_point(const real tval,const real pinf,
                        real & vpc_tpr, real & vpc_cav);
    real regula_falsi_kernel(real & tm, real & tp,
                             real & vpcm, real & vpcp,
                             real & vpc_tpr_m, real & vpc_tpr_p,
                             real & vpc_cav_m, real & vpc_cav_p,
                             const real pinf);

    Momentum & ns;
    EnthalpyFD & enthFD;
    PhaseChange4 & pc;
    CavityPressure & capr;
    AC & multigrid_cavity;
    TIF & tsat;
    const Times & time;

    Scalar & f;
    Scalar & press;
    Scalar & tpr;
    const Scalar & tpr_old;
    Vector & uvw;
    const Vector & uvw_old;
    Vector & uvw_cav;

    const std::function<real(const real t)> model;

    real t_new, t_old;

    Cycle c0, c1;
    ResRat rr_cav, rr_ns;
    ResTol rt_cav;
    std::array<MaxIter,3> mi;

    ResTol errtol;
};
#endif
