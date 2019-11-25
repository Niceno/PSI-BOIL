#ifndef MARCHING_SQUARES_H
#define MARCHING_SQUARES_H

#include "../heaviside.h"
#include "../../Topology/topology.h"

////////////////////////
//                    //
//  Marching Squares  //
//                    //
////////////////////////
class MarchingSquares : public Heaviside {
  public:
    MarchingSquares(const Comp & MCOMP, 
                    const Scalar * CLR, Scalar * PHI = NULL, 
                    Scalar * ADENS = NULL, const real CLRSURF = 0.5);
    ~MarchingSquares() {};

    virtual void evaluate_nodes();
    /* vf = area fraction in 2D terms */
    virtual real vf(const int i, const int j, const int k);
    /* ad = length density in 2D terms */
    virtual real ad(const int i, const int j, const int k);

    void topology(Topology & topo);

  protected:
    int construct_grid(const int i, const int j, const int k,
                       CELL2D & grid);

    real ad(const int i, const int j, const int k,
            std::vector<LINE> & lines);

    real line_density(const std::vector<LINE> & lines,
                      const real & surf);
    int extract_line_parameters(const std::vector<LINE> & lines,
                                real & nx, real & ny, real & nalpha);

    Scalar stmp;
    const Comp perpendicular;
    int ofx, ofy, ofz;
};

#endif
