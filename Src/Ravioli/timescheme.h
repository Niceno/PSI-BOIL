#ifndef TIMESCHEME_H
#define TIMESCHEME_H

#include "../Global/global_precision.h"

//////////////////
//              //
//  TimeScheme  //
//              //
//////////////////
class TimeScheme {

  public:
    TimeScheme(const real & a, const real & b, const real & c) 
     {N_ = a; N_m_1=b; N_m_2=c;}

    static const TimeScheme backward_euler() 
      {return TimeScheme(1.0, 0.0, 0.0);}

    static const TimeScheme forward_euler() 
      {return TimeScheme(0.0, 1.0, 0.0);}

    static const TimeScheme crank_nicolson() 
      {return TimeScheme(0.5, 0.5, 0.0);}

    static const TimeScheme adams_bashforth() 
      {return TimeScheme(0.0, 1.5,-0.5);}

    static const TimeScheme steady() 
      {return TimeScheme(1.0, 0.0, 0.0);}

    real N  () {return N_;}
    real Nm1() {return N_m_1;}
    real Nm2() {return N_m_2;}

  private:
    TimeScheme() {}
    real N_;
    real N_m_1;
    real N_m_2;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: timescheme.h,v 1.6 2010/03/25 08:15:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
