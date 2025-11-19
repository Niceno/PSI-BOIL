#ifndef LINEAR_H
#define LINEAR_H

#include "../../Parallel/mpi_macros.h"
#include <iostream>
#include <cmath>
#include "../../Ravioli/restol.h"
#include "../../Ravioli/resrat.h"
#include "../../Ravioli/maxiter.h"
#include "../../Domain/domain.h"
#include "../../Matrix/matrix.h"
#include "../../Field/Scalar/scalar.h"
#include "../Preconditioner/preconditioner.h"

/***************************************************************************//**
*  \brief Parent class for all linear algebraic solvers.
*
*  It stores all the scalar arrays (Scalar) needed for for solver algorithm,
*  as well as preconditioning matrix (Matrix) and Preconditioner. 
*  It also contains virtual function "solve", which calls solving algorithm 
*  implementation of a child solver.
*******************************************************************************/

//////////////
//          //
//  Linear  //
//          //
//////////////
class Linear {
  public:
    Linear(const Domain & d, const Prec & pc = Prec::di()); 

    //! Pointer to preconditioner associated with the solver. 
    Preconditioner * prec; 

    //! Calls solving algorythm implemented in the child.
    /*!
        \param A - system matrix,
        \param x - unknown vector,
        \param b - right hand side vector,
        \param iter - maximum number of iterations.
        \param var_name - variable name, used for timing routines (Timer)
        \param rt - residual tolerance, 
        \param rr - residual ratio.     
    */
    virtual bool solve(Matrix & A,              
                       Scalar & x,               
                       Scalar & b,                
                       const MinIter & mini,
                       const MaxIter & mi,
                       const char * var_name,
                       const ResRat & rr = ResRat(), 
                       const ResTol & rt = ResTol(),
                       const real scale = 1.0,
                       const int stalecount = -1,
                       const bool precform = true) = 0;

    //! Pointer to the domain on which the solver is defined.
    const Domain * domain() const {return dom;}

  protected:
    //! Array for child solver(s).
    Scalar p,  // cg, cgs, bicgs
           q,  // cg, cgs
           r,  // cg, cgs, bicgs
           s,  //     cgs, bicgs
           t,  //          bicgs
           u,  //     cgs     
           v,  //          bicgs
           z,  // cg     
           p_, //     cgs, bicgs
           q_, //     cgs
           r_, //     cgs, bicgs
           s_, //          bicgs
           u_, //     cgs
           v_; //     cgs

  private:
    const Domain * dom;
};

#endif
