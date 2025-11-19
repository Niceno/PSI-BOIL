#ifndef MAXITER_H
#define MAXITER_H

/***************************************************************************//**
*  \brief Ravioli class for safe parameter passing to linear solvers.
*
*  This class, when used as a parameter to linear solver, holds the value 
*  of maximum number of iterations for a solver. If, during a linear solver 
*  number of iterations exceeds this value, iterations are stopped. 
*******************************************************************************/

///////////////
//           //
//  MaxIter  // 
//           //
///////////////
class MaxIter {
  public:
    explicit MaxIter()              : val(20) {};
    explicit MaxIter(const int & v) : val(v)  {};
    operator int () const {return val;}

  private:
    int val;
};

#endif
