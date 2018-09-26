#ifndef POSITION_H
#define POSITION_H

#include "../../Ravioli/comp.h"

////////////////
//            //
//  Position  //
//            //
////////////////
class Position {
  public:
    Position()                          {xp[0] = 0;  xp[1] = 0;  xp[2] = 0;}
    Position(real xa, real ya, real za) {xp[0] = xa; xp[1] = ya; xp[2] = za;}

    void set_to(real xa, real ya, real za) {xp[0]=xa;   xp[1]=ya;   xp[2]=za;}
    void set_to(const Position & o) {for(int i=0; i<DIM; i++) xp[i] = o.xp[i];}

    real x() const {return xp[0];}
    real y() const {return xp[1];}
    real z() const {return xp[2];}

    real & xyz(const Comp & m) {return *(xp+~m);}

  private:
    real xp[DIM];
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: position.h,v 1.3 2013/08/12 05:39:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
