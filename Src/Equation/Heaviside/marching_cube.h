#ifndef MARCHING_CUBE_H
#define MARCHING_CUBE_H

#include "../../Field/Scalar/scalar.h"

/////////////////////
//                 //
//  Marching Cube  //
//                 //
/////////////////////
class MarchingCube {
  public:
    MarchingCube(const Scalar * CLR, const real CLRSURF = 0.5) {
      clr = CLR;
      clrsurf = CLRSURF;	
    }  
    ~MarchingCube() {}

    real volume(const int i, const int j, const int k);
    real area(const int i, const int j, const int k);
    real surface(const int dir, const Comp & mcomp,
                 const int i, const int j, const int k);

  protected:
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
    real clrsurf;
    
    real surfval(CELLFACE grid, real isolevel);
    FACEVAL faceval(const int dir, const Comp & mcomp,
                    const int i, const int j, const int k);

    XY VertexInterp2D(real isolevel, XY p1, XY p2, real valp1, real valp2);
    real SurfaceArea3(XY v1, XY v2, XY v3);
    real SurfaceArea4(XY v1, XY v2, XY v3, XY v4);
    real SurfaceArea5(XY v1, XY v2, XY v3, XY v4, XY v5);

};

#endif
