#ifndef ENTHALPYFDHIGHPR_H
#define ENTHALPYFDHIGHPR_H

#include "../enthalpyfd.h"

///////////////////////////////
//  for high Prandtl number  //
//  modify convection term   //
///////////////////////////////

class EnthalpyFDhighPr : public EnthalpyFD {
  public:
    EnthalpyFDhighPr(const Scalar & phi,
                     const Scalar & f,
                     const Scalar & clr,
                     const Vector & u,
                     Times & t,
                     Linear * sm,
                     Matter * flu,
                     const real tsat,
                     Matter * sol = NULL) :
    EnthalpyFD(phi,f,clr,u,t,sm,flu,tsat,sol) {}
    ~EnthalpyFDhighPr() {};
    void new_time_step(const Scalar * diff_eddy = NULL);
    void convection();
    void convection(Scalar * sca);

  protected:
    void setflag();
};
#endif
