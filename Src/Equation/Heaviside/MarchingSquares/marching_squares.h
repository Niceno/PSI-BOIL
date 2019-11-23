#ifndef MARCHING_SQUARES_H
#define MARCHING_SQUARES_H

#include "../heaviside.h"

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

  protected:
    int construct_grid(const int i, const int j, const int k,
                       CELL2D & grid);

    real line_density(const std::vector<LINE> & lines,
                      const real & surf);

    const Comp perpendicular;
    int ofx, ofy, ofz;
};

#endif
