#ifndef KRYLOV_H
#define KRYLOV_H

#include "../linear.h"

/***************************************************************************//**
*  \brief Parent class for all Krylov subspace solvers.
*******************************************************************************/

//////////////
//          //
//  Krylov  //
//          //
//////////////
class Krylov : public Linear {
  public:
    Krylov(const Domain & d, const Prec & pc = Prec::di())
      : Linear(d,pc) {}; 

    virtual bool solve(Matrix & A,              
                       Scalar & x,               
                       Scalar & b,                
                       const MaxIter & mi,
                       const char * var_name = NULL,
                       const ResRat & rr = ResRat(), 
                       const ResTol & rt = ResTol(),
                       const real scale = 1.0,
                       const int stalecount = -1,
                       const bool precform = true) = 0;
};

#endif

#include "cg.h"
#include "cgs.h"
#include "bicgs.h"
