#ifndef RESTOL_H
#define RESTOL_H

#include "../Global/global_precision.h"
#include "../Global/global_constants.h"

/***************************************************************************//**
*  \brief Ravioli class for safe parameter passing to linear solvers.
*
*  This class, when used as a parameter to linear solver, holds the value 
*  of residual tolerance for a solver. If, during a linear solver iterations 
*  residuals fall bellow this value, iterations are stopped. 
*
*  \note This criterion is probably more usefull for pressure correction 
*        equation.
*******************************************************************************/

//////////////
//          //
//  ResTol  // -> residual tolerance
//          //
//////////////
class ResTol {
  public:
    explicit ResTol()               : val(boil::femto) {};
    explicit ResTol(const real & v) : val(v)       {};
    operator real () const {return val;}

  private:
    real val;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: restol.h,v 1.4 2016/04/06 16:36:41 sato Exp $'/
+-----------------------------------------------------------------------------*/
