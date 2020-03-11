/***************************************************************************//**
*  \brief Direct solver by Gaussian elimination.
*
*  Direct solution of system of linear algebraic equations by Gaussian 
*  elimination. It was designed for coarsest levels of AdditiveCorrection 
*  CURRENTLY OBSOLETE!!!!
*******************************************************************************/

#ifndef GAUSS_H
#define GAUSS_H

#include "../../../Parallel/mpi_macros.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "../../../Global/global_malloc.h"
#include "../../../Matrix/matrix.h"

/////////////
//         //
//  Gauss  //
//         //
/////////////
class Gauss {

  public:
    Gauss() {}

    //! Solves the system of linear algebraic equations.
    /*!
        \param A - system matrix,
        \param X - unknow vector (variable being solved),
        \param B - right hand side vector. 
    */
    void solve(Matrix & A, Scalar & X, Scalar & B);

  private:
    void legs(real ** a, real * b, real * x, int * indx, int n);
    void elgs(real ** a, int * indx, int n);
};

#endif
