#ifndef ENTHALPYFD_H
#define ENTHALPYFD_H

#include "../../../Parallel/mpi_macros.h"
#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Solver/Gauss/gauss.h"
#include "../../../Timer/timer.h"

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
        \param f   - extarnal source array (\f$\dot{q}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Krylov subspace solver. It acts as a solver, or as a
                     smoother for AC multirid.
        \param flu - Holds all fluid properties (\f$\rho, C_p, \lambda\f$),
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
    */
    EnthalpyFD(const Scalar & phi, 
               const Scalar & f,
               const Scalar & clr,
               const Vector & u,
               Times & t,
               Krylov * sm,
               Matter * flu,
               const real tsat,
               Matter * sol = NULL);
    ~EnthalpyFD();

    void new_time_step(const Scalar * diff_eddy = NULL);
    void solve(const ResRat & fact, const char * name = NULL);
    void solve_sor(const int & it, const real & r, const char * name = NULL);
	  
    //! Direct solver introduced just for checking it.
    void direct() {
      for_ijk(i,j,k)
        fnew[i][j][k] = fold[i][j][k]
                      + cnew[i][j][k] * conv_ts.N()
                      + fbnd[i][j][k]
                      + fext[i][j][k];

      boil::timer.start("enthalpy solver");
      Gauss gs;
      gs.solve(A, phi, fnew);
      boil::timer.stop("enthalpy solver");
    }

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

    real get_turbP(){return turbP;}
    void set_turbP(real a){
      turbP=a;
      boil::oout<<"EnthalpyFD:turbP= "<<turbP<<"\n";
    }

    void convection();
  protected:
    void create_system(const Scalar * diff_eddy = NULL);
    void create_system_innertial();
    void create_system_diffusive(const Scalar * diff_eddy = NULL);
    void create_system_bnd();
    real update_rhs();
    void convection(Scalar * sca);
    void diffusion_fd(const Scalar * diff_eddy = NULL);
    real dVFD(const int i, const int j, const int k){
      real vol = 0.5 * (phi.dxw(i)+phi.dxe(i))
               * 0.5 * (phi.dys(j)+phi.dyn(j))
               * 0.5 * (phi.dzb(k)+phi.dzt(k));
      return vol;
    };
    void setflag();
    void diff_matrix(real & am, real & ac, real & ap
                , real & tm, real & tc, real & tp
                , real & aflagm, real & aflagp
                , const real vol, const real area
                , const bool onm, const bool onc, const bool onp
                , const bool ofm, const bool ofc, const bool ofp
                , const real lsm, const real lsc, const real lsp
                , const real lfm, const real lfc, const real lfp
                , const real clm, const real clc, const real clp
                , real dxm, real dxp
                , real fdm, real fdp, real fdms, real fdps
                , real pm, real pc, real pp
                , const real edm, const real edc, const real edp
                , const int i, const int j, const int k, const Comp m);

    const Scalar * clr;
    real tsat,rhol,rhov,cpl,cpv,lambdal,lambdav,clrsurf,epsl;
    bool store_clrold;
    Scalar clrold;
    ScalarInt iflag;
    real turbP; //turbulent Prandtl number
    bool laminar;

};
#endif
