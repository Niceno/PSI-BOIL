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
        \param sm  - Krylov subspace solver. It acts as a solver, or as a
                     smoother for AC multirid.
        \param mat - Holds fluid properties (\f$\rho\f$).
    */
    Pressure(const Scalar & phi, 
             const Scalar & f,  
             const Vector & v, 
             Times & t, 
             Krylov * sm,
             Matter * mat);
    ~Pressure();
	  
    //! Discretize the system of equations. 
    void discretize(); 

    //! Computes right hand side (velocity diverence) for pressure equation.
    real update_rhs();
 
    // ghost fluid method
    void ghost(const Scalar & c, const Scalar & k);
  protected:
};	

#endif
