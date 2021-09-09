#ifndef ENTHALPYFD_H
#define ENTHALPYFD_H

#include <cmath>
#include "../../../Parallel/mpi_macros.h"
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../CommonHeatTransfer/commonheattransfer.h"

/***************************************************************************//**
*  \brief Discretizes and solves enthalpy conservaion equation.
*
*  The eqation is discretized in integral form:                       
*  \f[
*        \int_V \frac{\partial \rho C_p T}{\partial t} dV
*      + \int_S \rho C_p {\bf u} T \, dS
*      = \int_S \lambda \nabla T \, dS
*      + \dot{Q}
*      \; \; \; \;
*      [\frac{J}{s} = W]
*  \f] 
*  where \f$T \; [K]\f$ is temperature, \f$\rho \; [\frac{kg}{m^3}]\f$ is 
*  density, \f$C_p \; [\frac{J}{kgK}]\f$ is thermal capacity, \f${t} \; [s]\f$ 
*  is time, \f${\bf u} \; [\frac{m}{s}]\f$ is convective velocity,
*  \f$\lambda \; [\frac{W}{mK}]\f$ is thermal conductivity and 
*  \f$\dot{Q} \; [\frac{J}{s}]\f$ is (external) heat source rate. 
*******************************************************************************/

////////////////
//            //
//  Enthalpy  //
//            //
////////////////
class EnthalpyFD : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - temperature (\f$T\f$),
        \param f   - external source array (\f$\dot{q}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Linear solver. It acts as a solver, or as a
                     smoother for AC multigrid.
        \param flu - Holds all fluid properties (\f$\rho, C_p, \lambda\f$),
        \param topo - properties of the free surface (from ITM).
        \param tifmodel - interfacial temperature model.
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
    */

    /* Most general constructor */
    EnthalpyFD(const Scalar & phi, 
               const Scalar & f,
               const Vector & umixed,
               const Vector & uliq,
               const Vector & ugas,
               Times & t,
               Linear * sm,
               Matter * flu,
               const CommonHeatTransfer & cht,
               Matter * sol = NULL);

    /* Delegating constructor */
    EnthalpyFD(const Scalar & phi,
               const Scalar & f,
               const Vector & umixed,
               Times & t,
               Linear * sm,
               Matter * flu,
               const CommonHeatTransfer & cht,
               Matter * sol = NULL) :
    EnthalpyFD(phi,f,umixed,umixed,umixed,t,sm,flu,
               cht,sol) {};

    void new_time_step(const Scalar * diff_eddy = NULL);
    void inertial(Scalar & sca, const bool interface_crossed, const Old old); 
    void inertial(const bool interface_crossed, const Old old); 
    void convective_time_step(Scalar & sca);
    void convective_time_step();
    virtual void convection();
    virtual void diffusion(const Scalar * diff_eddy = NULL) {
      evaluate_diffusion(Old::yes,diff_eddy);
      return;
    }
    virtual void solve(const ResTol & toler, const ResRat & fact,
                       const char * name = NULL);
    virtual void solve(const ResTol & toler, const char * name = NULL) {
      solve(toler,ResRat(-1.),name);
    } 
    virtual void solve(const ResRat & fact, const char * name = NULL) {
      solve(ResTol(boil::atto),fact,name);
    }

    //! Interface call to parent's discretization.
    void discretize(const Scalar * diff_eddy = NULL) {
      boil::timer.start("enthalpy discretize");
      create_system(diff_eddy);
      if(diff_eddy != NULL)laminar=false;
      boil::timer.stop("enthalpy discretize");
    }

#include "enthalpyfd_inline.h"

    /* unit tests */
    bool test_extrapolation(const int count);
    bool test_extrapolation(std::vector<real> & stencil,
                            const std::vector<real> & coefficients,
                            const std::vector<int> & cutpoints);

  protected:
    typedef real (EnthalpyFD::*coef_gen)(const real,const real,const real);

    void evaluate_diffusion(const Old old, const Scalar * diff_eddy = NULL);
    
    void create_system(const Scalar * diff_eddy = NULL);
    void create_system_innertial();
    void create_system_diffusive(const Scalar * diff_eddy = NULL) {
      evaluate_diffusion(Old::no,diff_eddy);
      return;
    }
    void create_system_bnd();
    real update_rhs();
    void convection(Scalar * sca);

    void cell_diffusion_fluid(const Comp m,
                          const int i, const int j, const int k,
                          const int ox, const int oy, const int oz,
                          const real x0,
                          const coef_gen coef_m, const coef_gen coef_p,
                          const ResistEval re, const Old old,
                          const real tscn, const real tscm,
                          std::vector<StencilPoint> & stencil,
                          real & Aw, real & Ac, real & Ae, real & F,
                          const Scalar * diff_eddy);

    void cell_diffusion_solid(const Comp m,
                          const int i, const int j, const int k,
                          const int ox, const int oy, const int oz,
                          const real x0,
                          const coef_gen coef_m, const coef_gen coef_p,
                          const real dSm, const real dSp,
                          const ResistEval re, const Old old,
                          const real tscn, const real tscm,
                          real & Aw, real & Ac, real & Ae, real & F,
                          const Scalar * diff_eddy);

    void kernel_fluid1(const std::array<ConnectType,3> & ctype,
                       const real cxm, const real cxp,
                       std::vector<StencilPoint> & stencil,
                       const real resinvm, const real resinvp,
                       real & Am, real & Ac, real & Ap, real & F);

    void kernel_fluid2(const std::array<ConnectType,3> & ctype,
                       const real cxm, const real cxp,
                       std::vector<StencilPoint> & stencil,
                       const real resinvm, const real resinvp,
                       const real reswallm, const real reswallp,
                       const std::array<real,3> resistvals,
                       const real dwsrcm, const real dwsrcp,
                       real & Am, real & Ac, real & Ap, real & F);

    void kernel_solid(const std::array<ConnectType,3> & ctype,
                      const real cxm, const real cxp,
                      const real pm, const real pp,
                      const std::array<real,3> resistvals,
                      const real dwsrcm, const real dwsrcp,
                      real & Am, real & Ac, real & Ap, real & F);

    /* needed in diffusion kernels */
    inline real resistance_multiplier(const real res0,
                                      const real res1) const {
      return res0 / (res0 + res1);
    }

    virtual real coef_x_m(const real dxm, const real dxp, const real x0);
    virtual real coef_x_p(const real dxm, const real dxp, const real x0);
    virtual real coef_y_m(const real dxm, const real dxp, const real x0);
    virtual real coef_y_p(const real dxm, const real dxp, const real x0);
    virtual real coef_z_m(const real dxm, const real dxp, const real x0);
    virtual real coef_z_p(const real dxm, const real dxp, const real x0);

    real face_value(const Sign matter_sig, const Comp m, const real vel,
                    const int i, const int j, const int k,
                    const int ofx, const int ofy, const int ofz,
                    const Old old);
    real extrapolate_value(const Sign & dir,
                           const std::vector<StencilPoint> & stencil,
                           const StencilPoint & ctm,
                           const StencilPoint & ctp,
                           const real xpos);
    void extrapolate_values(std::vector<StencilPoint> & stencil,
                            const StencilPoint & ctm, const StencilPoint & ctp);
    void point_extrapolation(std::vector<StencilPoint> & stencil,
                             const int i0, const int i1,
                             const std::vector<StencilPoint> & extcil);

    /* This points to solid if solid() = true and fluid() otherwise.
       So you can always dereference it without segfaults */
    const Matter * safe_solid; 

    const Vector * uliq, * ugas;
    Vector flux_liq, flux_gas;
    
    const CommonHeatTransfer & cht;
    /* gradt in convection difference order */
    AccuracyOrder ao_conv;

    const BndFlag bflag_struct;
    Scalar ftif;
    ScalarInt iflag,iflagold;
    bool laminar;

    /* reference ctypes */
    std::array<ConnectType,3> c_fff, c_sss,
                              c_iff, c_ffi, c_ifi,
                              c_sff, c_ffs, c_sfs,
                              c_fss, c_ssf, c_fsf,
                              c_iss, c_ssi, c_isi,
                              c_fsi, c_isf,
                              c_sfi, c_ifs;
};	
#endif
