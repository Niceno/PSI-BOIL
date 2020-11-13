#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "../linear.h"

/***************************************************************************//**
*  \brief Parent class for all iterative solvers. Currently, no preconditioning.
*******************************************************************************/

/////////////////
//             //
//  Iterative  //
//             //
/////////////////
class Iterative : public Linear {
  public:
    Iterative(const Domain & d)
      : Linear(d) {}; 

    virtual void solve(Matrix & A,              
                       Scalar & x,               
                       Scalar & b,                
                       const MaxIter & mi,
                       const char * var_name,
                       const ResRat & rr = ResRat(), 
                       const ResTol & rt = ResTol(),
                       const real scale = 1.0,
                       const int stalecount = -1,
                       const bool precform = true) = 0;
};

#endif

#include "gaussseidel.h"
#include "jacobi.h"
