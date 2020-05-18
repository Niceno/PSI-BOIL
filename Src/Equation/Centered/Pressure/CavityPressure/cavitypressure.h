#ifndef CAVITYPRESSURE_H
#define CAVITYPRESSURE_H

#include <functional> /* std::function */
#include "../pressure.h"
#include "../../../Topology/topology.h"

/***************************************************************************//**
*  \brief Discretizes and solves pressure-Poisson equation in cavity approx.
*
*  Liquid is considered incompressible. Gas is compressible but not solved.
*  For the moment, only one cavity is assumed.
*  See PhD thesis of Leon Malan (2017, under Zaleski) for details.
*******************************************************************************/

//////////////////////
//                  //
//  CavityPressure  //
//                  //
//////////////////////
class CavityPressure : public Pressure {
  public:
    //! Global constructor.
    CavityPressure(const Scalar & phi, 
                   const Scalar & f,  
                   const Vector & v, 
                   Times & t, 
                   Linear * sm,
                   Matter * liq,
                   Topology & topo,
                   const real tens,
                   const Scalar & curv,
                   Sign sig = Sign::pos());
    ~CavityPressure();

    inline real get_cavity_pressure() const { return cavity_pressure; }
    inline void set_cavity_pressure(const real cp) {
      cavity_pressure = cp;
      boil::oout<<"CavityPressure::value= "<<cavity_pressure<<"\n";
    }

    /* Written like this to allow to call parent discretize kernel */

    //! Discretize the system of equations.
    virtual void discretize(const Scalar * diff_eddy = NULL);

    //! Computes right hand side (velocity diverence) for pressure equation.
    virtual real update_rhs();
	  
  protected:
    void diff_matrix(const int i, const int j, const int k);

    ScalarInt iflag;
    Vector fs;
    Scalar clr, kappa;
    const Sign sig; /* pos: liquid is phi=1 and vice versa */

    real tens;
    real cavity_pressure;

    /* test for gas cells: initialized as a lambda in constructor */
    real clrsurf;
    std::function<bool(const real)> in_gas;
};	

#endif
