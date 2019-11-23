#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include "../heaviside.h"
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>

/////////////////////
//                 //
//  Marching Cube  //
//                 //
/////////////////////
class MarchingCubes : public Heaviside {
  public:
    MarchingCubes(const Scalar * CLR, Scalar * PHI = NULL, 
                 Scalar * ADENS = NULL,const real CLRSURF = 0.5);
    ~MarchingCubes() {};

    virtual real vf(const int i, const int j, const int k);
    virtual real area(const int i, const int j, const int k);
    real surface(const Sign sig, const Comp & mcomp,
                 const int i, const int j, const int k);

  protected:
    virtual void evaluate_nodes();

    real cellface_af(const Sign sig, const Comp & mcomp,
                     const int i, const int j, const int k);
    VAL2D cellface_covered(const Sign sig, const Comp & mcomp,
                           const int i, const int j, const int k);

    int find_vertices(const int & cubeindex,const CELL3D & grid,
                      const real & isolevel,TRIANGLE * triangles);
    real polygonise_area(const CELL3D & grid, const real & isolevel);
    real polygonise_volume(const CELL3D & grid, const real & isolevel);

    const std::array<const int,256> edgeTable;
    std::vector< std::vector<int> > triTable;
};

#endif
