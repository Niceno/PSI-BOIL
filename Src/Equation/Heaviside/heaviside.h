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

    real & value(const int i, const int j, const int k) 
      { return phi[i][j][k]; }
    real volume(const int i, const int j, const int k) 
      { return mc.volume(i,j,k); }
    real area  (const int i, const int j, const int k) 
      { return mc.area(i,j,k); }
    real area(const real p1, const real p2, const real p3,
              const real p4, const real p5, const real p6,
              const real p7, const real p8,
              const int i, const int j, const int k)
      {  return mc.area(p1,p2,p3,p4,p5,p6,p7,p8,
            i, j, k); }

   
  protected:
    const Domain * dom; 
    Scalar phi;
    MarchingCube mc;
};

#endif
