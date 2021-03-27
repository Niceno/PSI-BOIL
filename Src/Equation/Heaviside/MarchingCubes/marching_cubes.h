#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include "../heaviside.h"
#include <array>

//////////////////////
//                  //
//  Marching Cubes  //
//                  //
//////////////////////
class MarchingCubes : public Heaviside {
  public:
    MarchingCubes(const Scalar * CLR, Scalar * PHI = NULL, 
                 Scalar * ADENS = NULL, const real CLRSURF = 0.5);
    virtual ~MarchingCubes() {};
    //~MarchingCubes() {};

    virtual void evaluate_nodes();
    virtual int status(const int i, const int j, const int k);
    virtual real vf(const int i, const int j, const int k);
    virtual real ad(const int i, const int j, const int k);

    virtual void topology(Topology & topo, const real tol_wall, 
                          const bool use_interp, const bool use_subgrid) {
      boil::oout<<"MarchingCubes::topology: Underdevelopment! "
                <<"Exiting."<<boil::endl;
      exit(0);
    }

    virtual real surface(const Sign sig, const Comp & mcomp,
                         const int i, const int j, const int k);

    /* test approach for adens */
    virtual void marker_gradient(Scalar & mga, const bool evalflag = true);

  protected:
    int construct_grid(const int i, const int j, const int k,
                       CELL3D & grid);

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
