#ifndef ENTHALPYFDAXISYM_H
#define ENTHALPYFDAXISYM_H

#include "../enthalpyfd.h"

//////////////////
//              //
//  EnthalpyFD  //
//              //
//////////////////
class EnthalpyFDaxisym : public EnthalpyFD {
  public:
    EnthalpyFDaxisym(const Scalar & phi,
                     const Scalar & f,
                     const Vector & u,
                     const Vector & uliq,
                     const Vector & ugas,
                     Times & t,
                     Linear * sm,
                     Matter * flu,
                     const CommonHeatTransfer & cht,
                     Matter * sol = NULL) :

    EnthalpyFD(phi,f,u,uliq,ugas,t,sm,flu,cht,sol) {

      if(phi.domain()->is_cartesian()) {
           boil::oout<<"Warning: Initializing axisymmetric EnthalpyFD on a Cartesian "
                     <<"domain!"<<boil::endl;
         }
    }

    /* delegating constructor */
    EnthalpyFDaxisym(const Scalar & phi,
                     const Scalar & f,
                     const Vector & u,
                     Times & t,
                     Linear * sm,
                     Matter * flu,
                     const CommonHeatTransfer & cht,
                     Matter * sol = NULL) :
    EnthalpyFDaxisym(phi,f,u,u,u,t,sm,flu,cht,sol) {};

    ~EnthalpyFDaxisym() {};

  protected:

    virtual real coef_x_m(const real dxm, const real dxp, const real x0);
    virtual real coef_x_p(const real dxm, const real dxp, const real x0);
    virtual real coef_y_m(const real dxm, const real dxp, const real x0);
    virtual real coef_y_p(const real dxm, const real dxp, const real x0);
    virtual real coef_z_m(const real dxm, const real dxp, const real x0);
    virtual real coef_z_p(const real dxm, const real dxp, const real x0);

};
#endif
