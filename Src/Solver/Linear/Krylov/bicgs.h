#ifndef BICGS_H
#define BICGS_H

#include "krylov.h"

/***************************************************************************//**
*  \brief Krylov child embedding the Bi-Conjugate Gradient Stabilized (BiCGS) 
*         solver.
*******************************************************************************/

/////////////
//         //
//  BiCGS  //
//         //
/////////////
class BiCGS : public Krylov {
  public:
    BiCGS(const Domain & s, const Prec & pc) : Krylov(s, pc) {allocate(s);}
    BiCGS(const Domain & s)                  : Krylov(s)     {allocate(s);}

    virtual void solve(Matrix & A, Scalar & x, Scalar & b, 
                       const MaxIter & mi, const char * name = NULL,
                       const ResRat & rr = ResRat(),
                       const ResTol & rt = ResTol(),
                       const real scale = 1.0,
                       const int stalecount = -1,
                       const bool precform = true);

  private:
    void allocate(const Domain & d) {
      p. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      r .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      s. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      t .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      v .allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      p_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      r_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      s_.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }
};

#endif
