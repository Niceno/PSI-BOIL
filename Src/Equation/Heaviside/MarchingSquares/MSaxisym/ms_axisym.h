#ifndef MS_AXISYM_H
#define MS_AXISYM_H

#include "../marching_squares.h"

/////////////////////////////////////
//                                 //
//  Marching Axisymmetric Squares  //
//                                 //
/////////////////////////////////////
class MSaxisym : public MarchingSquares {
  public:
    MSaxisym(const Scalar * CLR, Scalar * PHI = NULL, 
             Scalar * ADENS = NULL, const real CLRSURF = 0.5) :
      MarchingSquares(Comp::j(),CLR,PHI,ADENS,CLRSURF) {};
    ~MSaxisym() {};

  protected:
    virtual real line_density(const std::vector<LINE> & lines,
                              const real & surf, const real & com);
};

#endif
