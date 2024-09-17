#ifndef ENTHALPYFDHIGHNU_H
#define ENTHALPYFDHIGHNU_H

#include "../enthalpyfd.h"

///////////////////////////////
//  for high Nusselt number  //
//  modify convection term   //
///////////////////////////////

class EnthalpyFDhighNu : public EnthalpyFD {
  public:
#if 1
    EnthalpyFDhighNu(const Scalar & phi,
                     const Scalar & f,
                     const Scalar & clr,
                     const Vector & u,
                     Times & t,
                     Linear * sm,
                     Matter * flu,
                     const real tsat,
                     Matter * sol = NULL) :
    EnthalpyFD(phi,f,clr,u,t,sm,flu,tsat,sol) {}
#endif
    ~EnthalpyFDhighNu() {};
    void new_time_step(const Scalar * diff_eddy = NULL);
    void convection();
    void convection(Scalar * sca);
    //void diffusion_fd(const Scalar * diff_eddy = NULL) {
    //}

  protected:
    void setflag();
  //private:
};

#endif
