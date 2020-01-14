#ifndef HEAVISIDE_H
#define HEAVISIDE_H

#include "../../Parallel/mpi_macros.h"
#include "../../Field/Scalar/scalar.h"
#include "../../Domain/domain.h"
#include "../../Global/global_realistic.h"
#include "../Topology/topology.h"
#include <vector>

/////////////////
//             //
//  Heaviside  //
//             //
/////////////////
class Heaviside { /* this class is an abstract class! */
  public:
    Heaviside(const Scalar * CLR, Scalar * PHI = NULL, Scalar * ADENS = NULL,
              const real CLRSURF = 0.5) : 
      clr(CLR), dom((*CLR).domain()), phi(PHI), adens(ADENS), clrsurf(CLRSURF) 
    {
      /* for vertex interpolation reasons */
      assert(boil::nano>boil::pico);
    }
    ~Heaviside() {};

    const Domain * domain() const {return dom;}
    void calculate(const bool evalflag = true);
    void calculate_vf(const bool evalflag = true);
    void calculate_adens(const bool evalflag = true);

    /* pure virtual functions */
    virtual void evaluate_nodes() = 0;
    virtual real ad(const int i, const int j, const int k) = 0;
    virtual real vf(const int i, const int j, const int k) = 0;

    virtual void topology(Topology & topo, const real tol_wall, 
                          const bool use_interp, const bool use_subgrid) = 0;
    void cal_fs_interp(const Scalar & scp,Vector & fs,
                       const real tol_wall, const bool use_subgrid);

    void fs_bnd_nosubgrid(const Scalar & scp, Vector & fs,
                          const real & tol_wall);
    void fs_bnd_geometric(const Scalar & scp, Vector & fs, 
                          const real & tol_wall);
    void fs_bnd_1D(const Scalar & scp, Vector & fs, 
                   const real & tol_wall, const Sign & sig);

    real operator() (const int i, const int j, const int k) const {
      return (*phi)[i][j][k];
    }

    Scalar nodalvals;

  protected:
    const Scalar * clr;
    const Domain * dom; 
    Scalar * phi;
    Scalar * adens;
    real clrsurf;

#include "heaviside_geometry.h"

    /* interpolation */
    XY  VertexInterp(const real & isolevel,
                     const XY & p1, const XY & p2, 
                     const real & valp1, const real & valp2);
    XYZ VertexInterp(const real & isolevel,
                     const XYZ & p1, const XYZ & p2,
                     const real & valp1, const real & valp2);

    /* cross product */
    XYZ CrossProduct(const XYZ & p1, const XYZ & p2);

    /* shoelace */
    real Shoelace(const XY & v1, const XY & v2, const XY & v3);
    real Shoelace(const XY & v1, const XY & v2, const XY & v3,
                  const XY & v4);
    real Shoelace(const XY & v1, const XY & v2, const XY & v3,
                  const XY & v4, const XY & v5);

    /* surface divergence */ 
    real triangle_surface_divergence(const TRIANGLE & t);

    /* main body of marching squares */
    real standing_square(const CELL2D & grid, const real & isolevel,
                         const real & totarea, std::vector<LINE> & lines);
  
    void cal_fs_geom(const Scalar & scp, 
                     const Scalar & nx, const Scalar & ny,
                     const Scalar & nz, const Scalar & nalpha,
                     Vector & fs,
                     const real tol_wall, const bool use_subgrid);
    
    real fs_val(const Comp m, const int i, const int j, const int k,
                const Scalar & nx, const Scalar & ny,
                const Scalar & nz, const Scalar & nalpha);
};

#endif
