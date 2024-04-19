#ifndef ENTHALPY_H
#define ENTHALPY_H

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
*        \rho C_p \int_V \frac{\partial T}{\partial t} dV
*      + \rho C_p ( \int_S {\bf u} T \, dS - T \int_S {\bf u} \, dS )
*      = \int_S \lambda \nabla T \, dS 
*      + \dot{Q}
*      \; \; \; \; 
*      [\frac{J}{s} = W]$
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
class Enthalpy : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - temperature (\f$T\f$),
        \param f   - extarnal source array (\f$\dot{q}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - linear solver. It acts as a solver, or as a
                     smoother for AC multigrid.
        \param flu - Holds all fluid properties (\f$\rho, C_p, \lambda\f$),
        \param sol - holds all solid properties (\f$\rho, C_p, \lambda\f$).
    */
    Enthalpy(const Scalar & phi, 
                const Scalar & f,
                const Vector & u, 
                Times & t,
                Linear * sm,
                Matter * flu,
                Matter * sol = NULL); 
    ~Enthalpy();
	  
    //! Direct solver introduced just for checking it.
    void direct() {
      for_ijk(i,j,k)
        fnew[i][j][k] = fold[i][j][k]; // + fext ....

      boil::timer.start("enthalpy solver");
      Gauss gs;
      gs.solve(A, phi, fnew);
      boil::timer.stop("enthalpy solver");
    }

    //! Discretization.
    void discretize(const Scalar * diff_eddy = NULL) {
      boil::timer.start("enthalpy discretize");
      create_system(diff_eddy);
      boil::timer.stop("enthalpy discretize");
    }

    //void new_time_step(const Scalar * diff_eddy);
    void new_time_step(const Scalar * diff_eddy = NULL);

    real get_turbPr(){return turbPr;}
    void set_turbPr(real a){
      turbPr=a;
      boil::oout<<"EnthalpyFD:turbPr= "<<turbPr<<"\n";
    }

    real turb_diff(const Scalar * diff_eddy, const int i, const int j, const int k);

    void set_fvm_convection(bool b){
      fvm_convection=b;
      if (fvm_convection) {
        boil::oout<<"enthalpy: use finite volume method for convection\n";
      } else {
        boil::oout<<"enthalpy: use 1st order finite difference method for convection\n";
      }
    }
    bool get_fvm_convection(){return fvm_convection;};


  protected:
    void create_system(const Scalar * diff_eddy = NULL);
    void create_system_innertial(const Property * f_prop,
                                 const Property * s_prop = NULL);
    void create_system_diffusive(const Property * f_prop,
                                 const Property * s_prop = NULL,
                                 const Scalar * diff_eddy = NULL);

    //void create_system_bnd(const Property * f_prop = NULL);

    //  evaluate_diffusion(Old::no,diff_eddy);
    //  return;
    //}
    //void create_system_bnd();
  
    void new_time_step(const Property * f_prop,
                       const Property * s_prop = NULL,
                       const Scalar * diff_eddy = NULL);

    // convection
    void convection(Scalar * conv, const Property * prop);

    real turbPr; // turbulent Prandtl number
    bool fvm_convection; // true: use FVM for convection, false: use FDM for convection
};	

#endif
