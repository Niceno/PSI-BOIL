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
                     const Vector & uliq,
                     const Vector & ugas,
                     Times & t,
                     Krylov * sm,
                     Matter * flu,
                     Topology & topo,
                     TIF & tifmodel,
                     Matter * sol = NULL) :
    EnthalpyFD(phi,f,clr,u,uliq,ugas,t,sm,flu,topo,tifmodel,sol) {

      if(phi.domain()->is_cartesian()) {
           boil::oout<<"Warning: Initializing axisymmetric EnthalpyFD on a Cartesian "
                     <<"domain!"<<boil::endl;
         }
    }

    /* delegating constructor */
    EnthalpyFDaxisym(const Scalar & phi,
                     const Scalar & f,
                     const Scalar & clr,
                     const Vector & u,
                     Times & t,
                     Krylov * sm,
                     Matter * flu,
                     Topology & topo,
                     TIF & tifmodel,
                     Matter * sol = NULL) :
    EnthalpyFDaxisym(phi,f,clr,u,u,u,t,sm,flu,topo,tifmodel,sol) {};

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
