#ifndef STAGGERED_H
#define STAGGERED_H

#include "../../Parallel/mpi_macros.h"
#include "../equation.h"
#include "../../Matter/matter.h"
#include "../../Timer/timer.h"

/////////////////
//             //
//  Staggered  //
//             //
/////////////////
class Staggered : public Equation {

  public:
    Staggered(const Domain * d, 
              const Vector & v,
              const Vector & f,
              Times & t,
              Matter * flu,
              Matter * sol,
              Krylov * sm) : Equation(d,&t,flu,sol,sm), u(&v), fext(&f),
                             fold(*d), fnew(*d), 
                             cold(*d), cnew(*d), gradp(*d) {
                set_indices();
                for_m(m) {
                  fext(m)  = v(m).shape();
                  fold(m)  = v(m).shape();
                  fnew(m)  = v(m).shape();
                  cold(m)  = v(m).shape();
                  cnew(m)  = v(m).shape();
                  gradp(m) = v(m).shape();
                }
              }

    //! Staggered cell volume.
    real dV(const Comp & m, const int i, const int j, const int k) const 
     {return u.dV(m,i,j,k);}

    const Vector & val() const {return u;}

  protected:
    int si(const Comp & m) const {return u.si(m);}
    int sj(const Comp & m) const {return u.sj(m);}
    int sk(const Comp & m) const {return u.sk(m);}
    int ei(const Comp & m) const {return u.ei(m);}
    int ej(const Comp & m) const {return u.ej(m);}
    int ek(const Comp & m) const {return u.ek(m);}

    /* connection dimensions needed for discretization.
       keep i,j,k here instead of x,y,z to be more general. */
    real dxw(const Comp & m, const int i) const {return u.dxw(m,i);}
    real dxe(const Comp & m, const int i) const {return u.dxe(m,i);}
    real dys(const Comp & m, const int j) const {return u.dys(m,j);}
    real dyn(const Comp & m, const int j) const {return u.dyn(m,j);}
    real dzb(const Comp & m, const int k) const {return u.dzb(m,k);}
    real dzt(const Comp & m, const int k) const {return u.dzt(m,k);}

    /* staggered cell-face areas needed for discretization */
    real dSx(const Comp & m, const int i, const int j, const int k) const
     {return u.dSx(m,i,j,k);}
    real dSy(const Comp & m, const int i, const int j, const int k) const
     {return u.dSy(m,i,j,k);}
    real dSz(const Comp & m, const int i, const int j, const int k) const
     {return u.dSz(m,i,j,k);}

    Vector u, fold, fnew, fext; 
    Vector    cold, cnew, gradp;

  private:
    void set_indices();
    void set_ranges(const Dir & d, const Comp & m, 
                    Range<int> * i, Range<int> * j, Range<int> * k);
};

#endif
