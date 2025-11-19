#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include "../../Domain/domain.h"
#include "../../Matrix/matrix.h"

////////////
//        //
//  Prec  // -> Ravioli type
//        //
////////////
class Prec {
  public:
    static const Prec no () {return Prec(-1);}
    static const Prec di () {return Prec( 0);}
    static const Prec ic0() {return Prec( 1);}
    static const Prec ic2() {return Prec( 2);}
    static const Prec ic3() {return Prec( 3);}

    operator int () const {return val;}
    
  private:
    /* prevent creation of new preconditioners */
    explicit Prec(const int i) {val = i;}
    int val;
};	

/***************************************************************************//**
*  \brief Parent class for all preconditioners.           
*
*  It holds only two pure virtual functions "form" and "solve" which call   
*  children's algorythms for forming of preconditioning matrices and 
*  solution algorythms.
*******************************************************************************/

//////////////////////
//                  //
//  Preconditioner  //
//                  //
//////////////////////
class Preconditioner {
  public:
    Preconditioner(const Domain & d); 

    //! Forms the preconditioning matrix.
    /*! \param A              - system matrix (from Linear solver). */
    virtual void form(const Matrix & A, const Scalar & x) = 0;

    //! Solves the preconditioning system.
    /*! \param z              - unknow vector,
        \param r              - right hand side vector. */
    virtual void solve(Scalar & z, const Scalar & r) = 0;

  protected:
    //! Preconditioning matrix.
    Scalar Mc;
    Scalar Mw;
    Scalar Ms;
    Scalar Mb;
    Scalar Mtw;
    Scalar Mts;
    Scalar Mnw;

    Preconditioner() {}
};

/***************************************************************************//**
*  \brief Diagonal Preconditioner.
*******************************************************************************/

//////////////////////////////////
//                              //
//  Diagonal <- Preconditioner  //
//                              //
//////////////////////////////////
class Diagonal : public Preconditioner {
  public:
    Diagonal(const Domain & d) {
      Mc.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }

    void form(const Matrix & A, const Scalar & x);
    void solve(Scalar & z, const Scalar & r);
};

/***************************************************************************//**
*  \brief Incomplete Cholesky Preconditioner with zero fill-in.
*******************************************************************************/

////////////////////////////////////////////
//                                        //
//  IncompleteCholesky <- Preconditioner  //
//                                        //
////////////////////////////////////////////
class IncompleteCholesky0 : public Preconditioner {
  public:
    IncompleteCholesky0(const Domain & d) {
      Mc.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mw.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Ms.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mb.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }

    void form(const Matrix & A, const Scalar & x);
    void solve(Scalar & z, const Scalar & r);
};

/***************************************************************************//**
*  \brief Incomplete Cholesky Preconditioner with two fill-ins.
*******************************************************************************/

////////////////////////////////////////////
//                                        //
//  IncompleteCholesky <- Preconditioner  //
//                                        //
////////////////////////////////////////////
class IncompleteCholesky2 : public Preconditioner {
  public:
    IncompleteCholesky2(const Domain & d) {
      Mc. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mw. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Ms. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mb. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mtw.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mts.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }

    void form(const Matrix & A, const Scalar & x);
    void solve(Scalar & z, const Scalar & r);
};

/***************************************************************************//**
*  \brief Incomplete Cholesky Preconditioner with three fill-ins.
*******************************************************************************/

////////////////////////////////////////////
//                                        //
//  IncompleteCholesky <- Preconditioner  //
//                                        //
////////////////////////////////////////////
class IncompleteCholesky3 : public Preconditioner {
  public:
    IncompleteCholesky3(const Domain & d) {
      Mc. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mw. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Ms. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mb. allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mtw.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mts.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
      Mnw.allocate(d.ni()+1, d.nj()+1, d.nk()+1);
    }

    void form(const Matrix & A, const Scalar & x);
    void solve(Scalar & z, const Scalar & r);
};

#endif
