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
                     const Scalar & clr,
                     const Vector & u,
                     Times & t,
                     Krylov * sm,
                     Matter * flu,
                     TIF & tifmodel,
                     Matter * sol = NULL,
                     const Vector * fs = NULL,
                     const Scalar * adens = NULL) :
    EnthalpyFD(phi,f,clr,u,t,sm,flu,tifmodel,sol,fs,adens) {

      if(phi.domain()->is_cartesian()) {
           boil::oout<<"Warning: Initializing axisymmetric EnthalpyFD on a Cartesian "
                     <<"domain!"<<boil::endl;
         }
    }

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