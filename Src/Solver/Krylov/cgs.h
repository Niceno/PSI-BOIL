#ifndef CGS_H
#define CGS_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include "../../Domain/domain.h"
#include "../../Matrix/matrix.h"
#include "krylov.h"

/***************************************************************************//**
*  \brief Krylov child embedding the Conjugate Gradient Squared (CGS) solver.
*******************************************************************************/

///////////
//       //
//  CGS  //
//       //
///////////
class CGS : public Krylov {
  public:
    CGS(const Domain & s, const Prec & pc) : Krylov(s, pc) {allocate(s);}
    CGS(const Domain & s)                  : Krylov(s)     {allocate(s);}

    void solve(Matrix & A, Scalar & x, Scalar & b, const MaxIter & mi, 
               const char * name = NULL,
               const ResRat & rr = ResRat(), const ResTol & rt = ResTol());

  private:
    void allocate(const Domain & d) {
      p. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      q. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      r .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      s .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      u .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      p_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      q_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      r_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      u_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      v_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: cgs.h,v 1.15 2014/08/06 07:44:24 sato Exp $'/
+-----------------------------------------------------------------------------*/
