#ifndef HEAVISIDE_H
#define HEAVISIDE_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Domain/domain.h"
#include "marching_cube.h"

/////////////////
//             //
//  Heaviside  //
//             //
/////////////////
class Heaviside {
  public:
    Heaviside(const Scalar & PHI, const Scalar * CLR,
              const real CLRSURF = 0.5) : 
               dom(PHI.domain()), phi(&PHI), mc(CLR, CLRSURF) {}  
    ~Heaviside() {};

    const Domain * domain() const {return dom;}
    void calculate();
   
  protected:
    const Domain * dom; 
    Scalar phi;
    MarchingCube mc;
};

#endif
