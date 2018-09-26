#ifndef MATRIX_H
#define MATRIX_H

#include "../Field/Scalar/scalar.h"

/***************************************************************************//**
*  \brief Class for storing linear system equation matrices.
*
*  Stores matrices for three-dimensional orthogonal grids using compass
*  notation (w-e, s-n, b-t). Matrices are stored with six diagonals 
*  corresponding to neighboring coefficients, main diagonal and inverse
*  of the diagonal. The class does not offer a lot of functionality, i.e.
*  constructor is the only routine which is really implemented and it
*  merely allocates the memory for the matrix. If it was not for memory
*  allocation, Matrix could have been a structure.
*******************************************************************************/

//////////////
//          //
//  Matrix  //
//          //
//////////////
class Matrix {
  public:
    Matrix(const Scalar & v) :
      c (v), 
      ci(v),
      w (v), 
      e (v), 
      s (v), 
      n (v), 
      b (v), 
      t (v) 
      {
        c  = v.shape();
        ci = v.shape();
        w  = v.shape();
        e  = v.shape();
        s  = v.shape();
        n  = v.shape();
        b  = v.shape();
        t  = v.shape();
      }	
    ~Matrix() {}
    
    Scalar c,  ///< Central coefficient.
           ci; ///< Inverted central coefficient.                
    Scalar w,  ///< West (i-1) coefficient.
           e,  ///< East (i+1) coefficient.
           s,  ///< South (j-1) coefficient.
           n,  ///< North (j+1) coefficient. 
           b,  ///< Bottom (k-1) coefficient.
           t;  ///< Top (k+1) coefficient.
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: matrix.h,v 1.14 2014/08/06 08:48:39 sato Exp $'/
+-----------------------------------------------------------------------------*/
