#ifndef SOR_H
#define SOR_H

#include "iterative.h"

/***************************************************************************//**
*  \brief Iterative child embedding the Successive Over-Relaxation (SOR) solver.
*******************************************************************************/

//////////////
//          //
//  SOR     //
//          //
//////////////
class SOR : public Iterative {
  public:
    SOR(const Domain & s) : Iterative(s)     {allocate(s);}

    void solve(Matrix & A, Scalar & x, Scalar & b, 
                       const MinIter & mini,
                       const MaxIter & mi, const char * name = NULL,
                       const ResRat & rr = ResRat(),
                       const ResTol & rt = ResTol());

  private:
    void allocate(const Domain & s) {
      r .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
    }
};
#endif
