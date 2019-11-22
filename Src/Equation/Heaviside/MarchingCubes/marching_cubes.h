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

    virtual real volume(const int i, const int j, const int k);
    virtual real area(const int i, const int j, const int k);
    real surface(const Sign sig, const Comp & mcomp,
                 const int i, const int j, const int k);

  protected:
    real surfval(CELL2D grid, real isolevel);
    VAL2D faceval(const Sign sig, const Comp & mcomp,
                    const int i, const int j, const int k);

    real SurfaceArea3(XY v1, XY v2, XY v3);
    real SurfaceArea4(XY v1, XY v2, XY v3, XY v4);
    real SurfaceArea5(XY v1, XY v2, XY v3, XY v4, XY v5);

    real polygonise_area(CELL3D grid, real isolevel);
    real polygonise_volume(CELL3D grid, real isolevel);

    real triangle_area(const TRIANGLE t);
    real triangle_vol_area(const TRIANGLE t);

    const std::array<const int,256> edgeTable_volume, edgeTable_area;
    std::vector< std::vector<int> > triTable;
};

#endif
