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

    real value(const Scalar & u);
    real value(const Vector & u, const Comp & m);

    real get_scalar(const Scalar & u);
    real get_vector(const Vector & u, const Comp & m);

    real get_grid(const Scalar & u, const Comp & m);
    real get_grid(const Vector & u, const Comp & m, const Comp & n);

  private:
    //Location() : dom(NULL), m_i(-1), m_j(-1), m_k(-1) {};
    Location() : dom(NULL), m_i(+10), m_j(-1), m_k(-1) {}; //boil::BW
   
    const char   * name;
    const Domain * dom;

    int m_i;
    int m_j;
    int m_k;

    int found_here;
};

#endif
