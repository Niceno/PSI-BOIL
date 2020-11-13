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
    virtual void diffusion(const Scalar * diff_eddy = NULL);
    virtual void solve(const ResTol & toler, const ResRat & fact,
                       const char * name = NULL);
    virtual void solve(const ResTol & toler, const char * name = NULL) {
      solve(toler,ResRat(-1.),name);
    } 
    virtual void solve(const ResRat & fact, const char * name = NULL) {
      solve(ResTol(-1.),fact,name);
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
    
    void create_system(const Scalar * diff_eddy = NULL);
    void create_system_innertial();
    void create_system_diffusive(const Scalar * diff_eddy = NULL);
    void create_system_bnd();
    real update_rhs();
    void convection(Scalar * sca);
    void diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , bool & aflagm, bool & aflagp
                , real & sourceterm
                , const real x0, const coef_gen coef_m, const coef_gen coef_p
                , const real vol, const real aream, const real areap
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const real lvm, const real lvc, const real lvp
                , const real llm, const real llc, const real llp
                , const int clm, const int clc, const int clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , const int i, const int j, const int k, const Comp m);

    inline real resistance_multiplier(const real res1, const real res2) const;

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

    /* faster diffusion matrix and ftif construction */
    bool accelerated_no_solid;

    const Vector * uliq, * ugas;
    Vector flux_liq, flux_gas;
    
    const CommonHeatTransfer & cht;
    /* gradt in convection difference order */
    AccuracyOrder ao_conv;

    const BndFlag bflag_struct;
    Scalar ftif;
    ScalarInt iflag,iflagold;
    bool laminar;
};	
#endif
