#ifndef LIMITER_H
#define LIMITER_H

#include "../Parallel/mpi_macros.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include "global_minmax.h"
#include "global_precision.h"
#include "../Parallel/Out/print.h"

/***************************************************************************//**
*  \brief Ravioli class for safer argument passing in Limiter.
*
*  Used for setting the right limiting function in Limiter. The names of
*  the member functions are self-explanatory.
*
*  \note It is closely connected with Limiter. If new schemes are added
*        to Limiter, ConvScheme has to be updated accordingly. 
*******************************************************************************/

//////////////////
//              //
//  ConvScheme  //
//              //
//////////////////
class ConvScheme {
  public:
    static const ConvScheme upwind  () {return ConvScheme(1);}
    static const ConvScheme central () {return ConvScheme(2);}
    static const ConvScheme minmod  () {return ConvScheme(3);}
    static const ConvScheme smart   () {return ConvScheme(4);}
    static const ConvScheme muscl   () {return ConvScheme(5);}
    static const ConvScheme superbee() {return ConvScheme(6);}
    static const ConvScheme van_leer() {return ConvScheme(7);}
    static const ConvScheme mc      () {return ConvScheme(8);}
    
    operator int () const {return val;}
    
    //! Prints the convection scheme name.
    friend std::ostream & 
      operator << (std::ostream & os, const ConvScheme & cs);

  private:
    /* prevent creation of new b.c. */
    explicit ConvScheme(const int i) {val = i;}
    int val;
};	

/***************************************************************************//**
*  \brief Class representing the convection scheme.             
*
*  Implements all convection schemes in the program. It can be used for 
*  Centered and Staggered (only Momemntum) data types. 
*
*  \note It is closely connected with ConvScheme. If new schemes are added
*        to Limiter, ConvScheme has to be updated accordingly. 
*******************************************************************************/

///////////////
//           //
//  Limiter  //
//           //
///////////////
class Limiter {
  public:
    //! Basic constructor. 
    Limiter(const ConvScheme & cs); 

    //! Sets the appropriate convection scheme.               
    void set(const ConvScheme & cs);

    //! Implementation of a limiter. 
    real limit(const real ucp,
               const real phim, const real phic, const real phip);

  private:
    real upwind  (const real r) const;
    real central (const real r) const;
    real minmod  (const real r) const;
    real smart   (const real r) const;
    real muscl   (const real r) const;
    real superbee(const real r) const;
    real van_leer(const real r) const;
    real mc      (const real r) const;

    /* pointer to limiting function */
    real (Limiter::*b)(const real r) const;
};

#endif
