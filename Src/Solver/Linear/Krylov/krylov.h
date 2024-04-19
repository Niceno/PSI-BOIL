#ifndef KRYLOV_H
#define KRYLOV_H

#include "../linear.h"
//#include "../../Parallel/mpi_macros.h"
//#include <iostream>
//#include <cmath>
//#include "../../Ravioli/restol.h"
//#include "../../Ravioli/resrat.h"
//#include "../../Ravioli/maxiter.h"
//#include "../../Domain/domain.h"
//#include "../../Matrix/matrix.h"
//#include "../../Field/Scalar/scalar.h"
//#include "../Preconditioner/preconditioner.h"

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

    virtual void solve(Matrix & A,              
                       Scalar & x,               
                       Scalar & b,                
                       const MinIter & minit,
                       const MaxIter & maxit,
                       const char * var_name = NULL,
                       const ResRat & rr = ResRat(), 
                       const ResTol & rt = ResTol()) = 0;

};

#endif

#include "cg.h"
#include "cgs.h"
#include "bicgs.h"
