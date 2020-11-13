#ifndef GS_H
#define GS_H

#include "iterative.h"

/***************************************************************************//**
*  \brief Iterative child embedding the Gauss-Seidel (GS) solver.
*******************************************************************************/

//////////////
//          //
//  Jacobi  //
//          //
//////////////
class GaussSeidel : public Iterative {
  public:
    GaussSeidel(const Domain & s) : Iterative(s)     {allocate(s);}

    virtual void solve(Matrix & A, Scalar & x, Scalar & b, 
                       const MaxIter & mi, const char * name = NULL,
                       const ResRat & rr = ResRat(),
                       const ResTol & rt = ResTol(),
                       const real scale = 1.0,
                       const int stalecount = -1,
                       const bool precform = true);

  private:
    void allocate(const Domain & s) {
      r .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
    }
};

#endif
