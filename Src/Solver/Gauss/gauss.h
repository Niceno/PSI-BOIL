/***************************************************************************//**
*  \brief Direct solver by Gaussian elimination.
*
*  Direct solution of system of linear algebraic equations by Gaussian 
*  elimination. It was designed for coarsest levels of AdditiveCorrection 
*  multigrid solver and therefore handles only relativelly small systems.
*
*  \note
*  May use class Board for plotting system matrices, introduced for debugging.
*  This calls to debugging plottings are disabled by commenting them out.
*******************************************************************************/

#ifndef GAUSS_H
#define GAUSS_H

#include "../../Parallel/mpi_macros.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "../../Global/global_malloc.h"
#include "../../Matrix/matrix.h"
#include "../../Board/board.h"

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
    void plot_system(real ** a, real * b, int c, int nj, int nk, char * name);
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: gauss.h,v 1.8 2014/08/06 07:43:13 sato Exp $'/
+-----------------------------------------------------------------------------*/
