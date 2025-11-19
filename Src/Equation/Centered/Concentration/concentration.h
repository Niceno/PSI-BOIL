#ifndef CONCENTRATION_H
#define CONCENTRATION_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"

/***************************************************************************//**
*  \brief Discretizes and solves species conservaion equation.
*
*  The eqation is discretized in integral form:                       
*  \f[
*        \int_V \frac{\partial \rho \alpha}{\partial t} dV
*      + \int_S \rho {\bf u} \alpha \, dS
*      = \int_S \gamma \nabla \alpha \, dS
*      + \dot{M}
*      \; \; \; \;
*      [\frac{kg}{s}]
*  \f] 
*  where \f$\alpha \; [1]\f$ is species concentration varying from 0 to 1,
*  \f$\rho \; [\frac{kg}{m^3}]\f$ is mass density, \f${t} \; [s]\f$ is time, 
*  \f${\bf u} \; [\frac{m}{s}]\f$ is convective velocity,
*  \f$\gamma \; [\frac{kg}{ms}]\f$ is species diffusivity and 
*  \f$\dot{M} \; [\frac{kg}{s}]\f$ is (external) mass source rate.
*******************************************************************************/

/////////////////////
//                 //
//  Concentration  //
//                 //
/////////////////////
class Concentration : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - species concentrationi array (\f$\alpha\f$),
        \param f   - extarnal source array (\f$\dot{m}\f$),
        \param u   - convection velocity (\f${\bf u}\f$),
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm - Krylov subspace solver. It acts as a solver, or as a
                    smoother for AC multirid.
        \param flu - Holds all material properties (\f$\rho, \gamma\f$).
    */
    Concentration(const Scalar & phi, 
                  const Scalar & f,
                  const Vector & u, 
                  Times & t,
                  Krylov * sm,
                  Matter * flu);
    ~Concentration();
	  
    //! Interface call to parent's discretization.
    void discretize(const Scalar * eddy = NULL);
};	

#endif
