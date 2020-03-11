#ifndef JACOBI_H
#define JACOBI_H

#include "iterative.h"

/***************************************************************************//**
*  \brief Iterative child embedding the Jacobi solver.
*******************************************************************************/

//////////////
//          //
//  Jacobi  //
//          //
//////////////
class Jacobi : public Iterative {
  public:
    Jacobi(const Domain & s) : Iterative(s)     {allocate(s);}

    virtual void solve(Matrix & A, Scalar & x, Scalar & b, 
                       const MaxIter & mi, const char * name = NULL,
                       const ResRat & rr = ResRat(),
                       const ResTol & rt = ResTol());

  private:
    void allocate(const Domain & s) {
      r .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
    }
};

#endif
