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
    Heaviside(const Scalar * CLR, Scalar * PHI = NULL, Scalar * ADENS = NULL,
              const real CLRSURF = 0.5) : 
               dom((*CLR).domain()), phi(PHI), adens(ADENS),
               mc(CLR, dom, CLRSURF) {}  
    ~Heaviside() {};

    const Domain * domain() const {return dom;}
    void calculate();
    void calculate_heaviside();
    void calculate_adens();

    real volume(const int i, const int j, const int k) 
      { return mc.volume(i,j,k); }
    real area  (const int i, const int j, const int k) 
      { return mc.area(i,j,k); }
#if 0
    real area(const real p1, const real p2, const real p3,
              const real p4, const real p5, const real p6,
              const real p7, const real p8,
              const int i, const int j, const int k)
      {  return mc.area(p1,p2,p3,p4,p5,p6,p7,p8,
            i, j, k); }
#endif
    real operator() (const int i, const int j, const int k) {
    //const real operator() (const int i, const int j, const int k) {
      return (*phi)[i][j][k];
    }

  protected:
    const Domain * dom; 
    Scalar * phi;
    Scalar * adens;
    MarchingCube mc;
};

#endif
