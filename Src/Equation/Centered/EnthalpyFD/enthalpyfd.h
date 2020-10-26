#ifndef ENTHALPYFD_H
#define ENTHALPYFD_H

#include "../../../Parallel/mpi_macros.h"
#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Timer/timer.h"
#include "../../../Global/global_realistic.h"
#include "../../Tifmodel/tif.h"
#include "../../Topology/topology.h"
#include "../../../Ravioli/htwallmodel.h"

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
               Topology * topo,
               TIF & tifmodel,
               Matter * sol = NULL,
               HTWallModel * htwallmodel = NULL);

    /* Delegating constructor */
    EnthalpyFD(const Scalar & phi,
               const Scalar & f,
               const Vector & umixed,
               Times & t,
               Linear * sm,
               Matter * flu,
               Topology * topo,
               TIF & tifmodel,
               Matter * sol = NULL,
               HTWallModel * htwallmodel = NULL) :
    EnthalpyFD(phi,f,umixed,umixed,umixed,t,sm,flu,
               topo,tifmodel,sol,htwallmodel) {};


    ~EnthalpyFD();

    void new_time_step(const Scalar * diff_eddy = NULL);
    void convective_time_step();
    void explicit_diffusion(const Scalar * diff_eddy = NULL);
    void solve(const ResRat & fact, const char * name = NULL);
    void solve_sor(const int & it, const real & r, const char * name = NULL);

    real hflux_wall(const Scalar & s, const Dir d
                  , const Scalar * diff_eddy = NULL);
    real hflux_wall_ib(const Scalar * diff_eddy = NULL);

    //! Interface call to parent's discretization.
    void discretize(const Scalar * diff_eddy = NULL) {
      boil::timer.start("enthalpy discretize");
      create_system(diff_eddy);
      if(diff_eddy != NULL)laminar=false;
      boil::timer.stop("enthalpy discretize");
    }

    inline real get_turbP() const { return turbP; }
    inline void set_turbP(const real a) {
      turbP = a;
      boil::oout<<"EnthalpyFD::turbP= "<<turbP<<"\n";
    }

    inline bool get_no_solid_acceleration() const
      { return accelerated_no_solid; }
    inline void set_no_solid_acceleration(const bool flag) {
      accelerated_no_solid = flag;
      boil::oout<<"EnthalpyFD::no_solid_acceleration= "
                <<accelerated_no_solid<<"\n";
    }

    inline HTWallModel & heat_transfer_wall_model() {
      return *htwallmodel;
    }

    void convection();

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
                , real & aflagm, real & aflagp
                , real & sourceterm
                , const real x0, const coef_gen coef_m, const coef_gen coef_p
                , const real vol, const real aream, const real areap
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const int clm, const int clc, const int clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m);

    inline real resistance_multiplier(const real dx1, const real dx2,
                                      const real l1, const real l2,
                                      const real resistplus = 0.0) const;

    virtual real coef_x_m(const real dxm, const real dxp, const real x0);
    virtual real coef_x_p(const real dxm, const real dxp, const real x0);
    virtual real coef_y_m(const real dxm, const real dxp, const real x0);
    virtual real coef_y_p(const real dxm, const real dxp, const real x0);
    virtual real coef_z_m(const real dxm, const real dxp, const real x0);
    virtual real coef_z_p(const real dxm, const real dxp, const real x0);

    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);
    bool interface_old(const Sign dir, const Comp m,
                       const int i, const int j, const int k);

    real Tint(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k);
    real Tint_old(const int dir, const Comp &mcomp, const real frac,
              const int i, const int j, const int k);
    real Tint(const int i, const int j, const int k);
    real Tint_old(const int i, const int j, const int k);

    /* cell distances */
    real distance_center(const Sign sig, const Comp & m,
                         const int i, const int j, const int k);
    real distance_face(const Sign sig, const Comp & m,
                       const int i, const int j, const int k);

    /* new distances */
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

    /* old distances */
    real distance_int_old(const Sign dir, const Comp & m,
                          const int i, const int j, const int k,
                          real & tint);
    real distance_int_x_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint);
    real distance_int_y_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint);
    real distance_int_z_old(const Sign dir,
                            const int i, const int j, const int k,
                            real & tint);

    /* ghost distance */
    real ghost_distance(const Comp & m, const Sign & cell_marker,
                        const int i, const int j, const int k);

    real gradt_ib(const int dir, const Comp & mcomp,
                  const int i, const int j, const int k);

    /* conductivity function */
    real lambda(const int i, const int j, const int k,
                const Scalar * diff_eddy = NULL) const;
    real lambda_inv(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const;

    /* This points to solid if solid() = true and fluid() otherwise.
       So you can always dereference it without segfaults */
    const Matter * safe_solid; 

    /* faster diffusion matrix and ftif construction */
    bool accelerated_no_solid;

    const Vector * uliq, * ugas;

    TIF & tifmodel;  
    Topology * topo;
    HTWallModel * htwallmodel;

    real rhol,rhov,cpl,cpv,lambdal,lambdav,epsl;
    Scalar ftif;
    ScalarInt iflag,iflagold;
    real turbP; /* turbulent Prandtl number */
    bool laminar;
    bool default_value_for_htwallmodel;
};	
#endif
