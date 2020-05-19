#ifndef CAVITYPRESSURE_H
#define CAVITYPRESSURE_H

#include <functional> /* std::function */
#include <cmath>
#include "../centered.h"
#include "../../Topology/topology.h"
#include "../../../Parallel/communicator.h"
#include "../../../Parallel/Out/out.h"

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
class CavityPressure : public Centered {
  public:
    //! Global constructor.
    CavityPressure(const Scalar & phi, 
                   const Scalar & f,  
                   const Vector & v, 
                   Times & t, 
                   Linear * sm,
                   Matter * liq,
                   Topology * topo,
                   const real tens,
                   const Scalar & curv,
                   Sign sig = Sign::pos());
    ~CavityPressure();

    inline real get_cavity_pressure() const { return cavity_pressure; }
    inline void set_cavity_pressure(const real cp) {
      cavity_pressure = cp;
      boil::oout<<"CavityPressure::value= "<<cavity_pressure<<"\n";
    }

    //! Discretize the system of equations.
    virtual void discretize(const Scalar * diff_eddy = NULL);

    //! Computes right hand side (velocity diverence) for pressure equation.
    virtual real update_rhs();
	  
  protected:
    void diff_matrix(const int i, const int j, const int k);
    virtual void create_system_bnd();

    bool interface(const Sign dir, const Comp m,
                   const int i, const int j, const int k);

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

    real Pint(const int i, const int j, const int k);

    virtual real coef_x_m(const real dxm, const real dxp, const real x0);
    virtual real coef_x_p(const real dxm, const real dxp, const real x0);
    virtual real coef_y_m(const real dxm, const real dxp, const real x0);
    virtual real coef_y_p(const real dxm, const real dxp, const real x0);
    virtual real coef_z_m(const real dxm, const real dxp, const real x0);
    virtual real coef_z_p(const real dxm, const real dxp, const real x0);

    Topology * topo;
    ScalarInt iflag;
    Vector fs;
    Scalar kappa;
    const Sign sig; /* pos: liquid is phi=1 and vice versa */

    real tens;
    real cavity_pressure;

    /* test for gas cells: initialized as a lambda in constructor */
    std::function<bool(const int,const int,const int)> in_gas;
};	

#endif
