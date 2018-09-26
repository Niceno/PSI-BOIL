#ifndef LOCATION_H
#define LOCATION_H

#include "../monitor.h"

/***************************************************************************//**
*  \brief Defines a monitoring location (a point) in computational domain.
* 
*  It is particularly usefull to follow history of solution at certain 
*  location, during and unsteady run. For steady runs, it might be usefull
*  to check if solution stopped changing at certain position. 
*******************************************************************************/

////////////////
//            //
//  Location  //
//            //
////////////////
class Location : public Monitor {
  public:
    Location(const char *,
             const Domain & dom, const int I, 
                                 const int J, 
                                 const int K);

    void print(const Scalar & u);
    void print(const Vector & u, const Comp & m);

  private:
    Location() : dom(NULL), m_i(-1), m_j(-1), m_k(-1) {};
   
    const char   * name;
    const Domain * dom;

    int m_i;
    int m_j;
    int m_k;

    int found_here;
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: location.h,v 1.9 2010/11/06 20:18:50 niceno Exp $'/
+-----------------------------------------------------------------------------*/
