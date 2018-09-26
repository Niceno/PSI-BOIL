#ifndef CG_H
#define CG_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include "../../Domain/domain.h"
#include "../../Matrix/matrix.h"
#include "krylov.h"

/***************************************************************************//**
*  \brief Krylov child embedding the Conjugate Gradient (CG) solver.
*******************************************************************************/

//////////
//      //
//  CG  //
//      //
//////////
class CG : public Krylov {
  public:
    CG(const Domain & s, const Prec & pc) : Krylov(s, pc) {allocate(s);}
    CG(const Domain & s)                  : Krylov(s)     {allocate(s);}

    void solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi, 
               const char * name = NULL,
               const ResRat & rr = ResRat(), const ResTol & rt = ResTol());

  private:
    void allocate(const Domain & s) {
      p .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
      q .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
      r .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
      z .allocate(s.ni()+1, s.nj()+1, s.nk()+1);
    }
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: cg.h,v 1.15 2014/08/06 07:44:24 sato Exp $'/
+-----------------------------------------------------------------------------*/
