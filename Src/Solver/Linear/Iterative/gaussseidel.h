#ifndef GS_H
#define GS_H

#include "iterative.h"

/***************************************************************************//**
*  \brief Iterative child embedding the Gauss-Seidel (GS) solver.
*******************************************************************************/

///////////////////
//               //
//  GaussSeidel  //
//               //
///////////////////
class GaussSeidel : public Iterative {
  public:
    GaussSeidel(const Domain & s) : Iterative(s)     {allocate(s);}

    void solve(Matrix & A, Scalar & x, Scalar & b, 
                       const MinIter & minit,
                       const MaxIter & maxit, const char * name = NULL,
                       const ResRat & rr = ResRat(),
                       const ResTol & rt = ResTol());

  private:
    void allocate(const Domain & s) {
      r .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
    }
};

#endif
