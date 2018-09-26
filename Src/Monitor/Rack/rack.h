#ifndef RACK_H
#define RACK_H

#include "../Location/location.h"

/***************************************************************************//**
*  \brief Defines a rack of monitoring locations in a computational domain.
* 
*  Defined for a rack of points, along "i", "j" or "k" direction.
*  It's main purpose is to print the profile of certain variable at the
*  end of the run. 
*******************************************************************************/

////////////
//        //
//  Rack  //
//        //
////////////
class Rack : public Monitor {
  public:
    Rack(const char *,
         const Domain & dom, const Range<int> ri, 
                             const int j, 
                             const int k);
    Rack(const char *,
         const Domain & dom, const int i, 
                             const Range<int> rj, 
                             const int k);
    Rack(const char *,
         const Domain & dom, const int i, 
                             const int j, 
                             const Range<int> rk);

    void print(const Scalar & u);
    void print(const Vector & u, const Comp & m);

  private:
    const char   * name;
    std::vector<Location *> mons; 
    Range<int> r_i;
    Range<int> r_j;
    Range<int> r_k;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: rack.h,v 1.7 2009/05/09 19:06:48 niceno Exp $'/
+-----------------------------------------------------------------------------*/
