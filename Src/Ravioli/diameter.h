#ifndef DIAMETER_H
#define DIAMETER_H

#include <vector>

////////////////
//            //
//  Diameter  //
//            //
////////////////
class Diameter {
  public:
    Diameter()       {di = 0;}
    Diameter(real d) {di = d;}

    real d() const {return di;}

  private:
    real di;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: diameter.h,v 1.1 2013/07/31 14:22:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
