#ifndef MARCHING_CUBE_H
#define MARCHING_CUBE_H

#include "../../Field/Scalar/scalar.h"
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>

/////////////////////
//                 //
//  Marching Cube  //
//                 //
/////////////////////
class MarchingCube {
  public:
    MarchingCube(const Scalar * CLR, const Domain * dom, const real CLRSURF = 0.5);
    ~MarchingCube() {}

    real volume(const int i, const int j, const int k);
    real area(const int i, const int j, const int k);
    real surface(const int dir, const Comp & mcomp,
                 const int i, const int j, const int k);
 
    real area(const real p1, const real p2, const real p3,
              const real p4, const real p5, const real p6,
              const real p7, const real p8,
              const int i, const int j, const int k);

  protected:
    const Domain * domain() const {return dom;}

    /* 3D */
    typedef struct {
      real x,y,z;
    } XYZ;
    typedef struct {
       XYZ p[3];
       XYZ v[3];
    } TRIANGLE;
    typedef struct {
       XYZ p[8];
       real val[8];
    } GRIDCELL;
    typedef struct {
      XYZ v;
      XYZ ref;
      real refval;
    } VERT;

    /* 2D */
    typedef struct {
      real x,y;
    } XY;
    typedef struct {
      XY p[4];
      real val[4];
    } CELLFACE;
    typedef struct {
      CELLFACE face;
      real area;
    } FACEVAL;

    const Scalar * clr;
    const Domain * dom; 
    real clrsurf;
    
    real surfval(CELLFACE grid, real isolevel);
    FACEVAL faceval(const int dir, const Comp & mcomp,
                    const int i, const int j, const int k);

    XY VertexInterp2D(real isolevel, XY p1, XY p2, real valp1, real valp2);
    real SurfaceArea3(XY v1, XY v2, XY v3);
    real SurfaceArea4(XY v1, XY v2, XY v3, XY v4);
    real SurfaceArea5(XY v1, XY v2, XY v3, XY v4, XY v5);

    XYZ VertexInterp3D(real isolevel, XYZ p1, XYZ p2, real valp1, real valp2);
    XYZ CrossProduct(const XYZ p1, const XYZ p2);
    real DotProduct(const XYZ p1, const XYZ p2);
    XYZ NegateXYZ(const XYZ p);
    XYZ PlusXYZ(const XYZ p1, const XYZ p2);

    real polygonise_area(GRIDCELL grid, real isolevel);
    real polygonise_volume(GRIDCELL grid, real isolevel);

    real triangle_area(const TRIANGLE t);
    real triangle_vol_area(const TRIANGLE t);

    const std::array<const int,256> edgeTable_volume, edgeTable_area;
#if 0
    const std::array< const std::array<const int,16> , 256> triTable;
#else
    std::vector< std::vector<int> > triTable;
#endif 
};

#endif
