#ifndef MONITOR_H
#define MONITOR_H

#include "../Domain/domain.h"
#include "../Field/Vector/vector.h"
#include "../Ravioli/comp.h"

/***************************************************************************//**
*  \brief Defines monitoring positions in computational domain
*
*  They can be used for monitoring values of dependent variables during,
*  or at the enf of solution procedure. 
*******************************************************************************/

///////////////
//           //
//  Monitor  //
//           //
///////////////
class Monitor {
  public:
    //! Prints the value of Scalar at specified position(s).
    virtual void print(const Scalar & u) = 0;

    //! Prints the value of a Vector component at specified position(s).
    /*! 
        \param u - Vector to be printed,
        \param m - specifies the Vector component.
    */
    virtual void print(const Vector & u, const Comp & m) = 0;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: monitor.h,v 1.6 2009/05/09 19:06:48 niceno Exp $'/
+-----------------------------------------------------------------------------*/
