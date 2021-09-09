#ifndef CPAXISYM_H
#define CPAXISYM_H

#include "../cavitypressure.h"

////////////////////////////////////
//                                //
//  Axisymmetric Cavity Pressure  //
//                                //
////////////////////////////////////
class CPaxisym : public CavityPressure {
  public:
    CPaxisym(const Scalar & phi,
             const Scalar & f,
             const Vector & v,
             Times & t,
             Linear * sm,
             Matter * liq,
             Topology * topo,
             const Property * tens = NULL,
             const Scalar * curv = NULL,
             Sign sig = Sign::pos()) :
    CavityPressure(phi,f,v,t,sm,liq,topo,tens,curv,sig) {

      if(phi.domain()->is_cartesian()) {
           boil::oout<<"Warning: Initializing axisymmetric Cavity Pressure "
                     <<"on a Cartesian domain!"<<boil::endl;
         }
    }

    ~CPaxisym() {};

  protected:

    virtual real coef_x_m(const real dxm, const real dxp, const real x0);
    virtual real coef_x_p(const real dxm, const real dxp, const real x0);
    virtual real coef_y_m(const real dxm, const real dxp, const real x0);
    virtual real coef_y_p(const real dxm, const real dxp, const real x0);
    virtual real coef_z_m(const real dxm, const real dxp, const real x0);
    virtual real coef_z_p(const real dxm, const real dxp, const real x0);

};
#endif
