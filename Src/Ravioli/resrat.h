#ifndef RESRAT_H
#define RESRAT_H

#include "../Global/global_precision.h"
#include "../Global/global_constants.h"

/***************************************************************************//**
*  \brief Ravioli class for safe parameter passing to linear solvers.
*
*  This class, when used as a parameter to linear solver, holds the value 
*  of residual ration for a solver. If, during a linear solver iterations, 
*  ratio between current and initial residuals fall bellow this value, 
*  iterations are stopped. 
*
*  \note This criterion is probably more usefull for transport equations (all 
*        except pressure).
*******************************************************************************/

//////////////
//          //
//  ResRat  // -> residual ratio
//          //
//////////////
class ResRat {
  public:
    explicit ResRat()               : val(0.01) {};
    explicit ResRat(const real & v) : val(v)    {};
    operator real () const {return val;}

  private:
    real val;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: resrat.h,v 1.1 2011/05/25 11:28:03 niceno Exp $'/
+-----------------------------------------------------------------------------*/
