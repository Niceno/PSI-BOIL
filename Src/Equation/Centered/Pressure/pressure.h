#ifndef PRESSURE_H
#define PRESSURE_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Parallel/Out/out.h"

/***************************************************************************//**
*  \brief Discretizes and solves pressure-Poisson equation.
*
*  The eqation is discretized in integral form:                       
*  \f[
*        \int_S \frac{\nabla p}{\rho} \, dS
*      =          
*        \frac{1}{\Delta t}
*        \int_S {\bf u} \, dS
*      \; \; \; \;
*      [\frac{m^3}{s^2}]
*  \f] 
*  where \f$p \; [\frac{kg}{m s^2}]\f$ is pressure,
*  \f$\rho \; [\frac{kg}{m^3}]\f$ is density, 
*  \f${\bf u} \; [\frac{m}{s}]\f$ is velocity, and
*  \f$\Delta {t} \; [s]\f$ is numerical time step.
*******************************************************************************/

////////////////
//            //
//  Pressure  //
//            //
////////////////
class Pressure : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - pressure correction (\f$p\f$),
        \param f   - right hand side array,
        \param v   - velocity field, 
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Linear solver. It acts as a solver, or as a
                     smoother for AC multigrid.
        \param mat - Holds fluid properties (\f$\rho\f$).
    */
    Pressure(const Scalar & phi, 
             const Scalar & f,  
             const Vector & v, 
             Times & t, 
             Linear * sm,
             Matter * mat);
    ~Pressure();

    /* Written like this to allow child to call parent discretize kernel */
	  
    //! Discretize the system of equations. 
    virtual void discretize(const Scalar * diff_eddy = NULL) {
      discretize_pressure(diff_eddy);
    }

    //! Computes right hand side (velocity diverence) for pressure equation.
    virtual real update_rhs() {
      return update_rhs_pressure();
    }
 
    // ghost fluid method
    void ghost(const Scalar & c, const Scalar & k);
  protected:
     
    void discretize_pressure(const Scalar * diff_eddy = NULL);
    real update_rhs_pressure();

};	

#endif
