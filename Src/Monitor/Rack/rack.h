#ifndef RACK_H
#define RACK_H

#include <memory>
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
    /* 1D */
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
    /* 2D */
    Rack(const char *,
         const Domain & dom, const int i,
                             const Range<int> rj,
                             const Range<int> rk);
    Rack(const char *,
         const Domain & dom, const Range<int> ri,
                             const int j,
                             const Range<int> rk);
    Rack(const char *,
         const Domain & dom, const Range<int> ri,
                             const Range<int> rj,
                             const int k);

    void print(const Scalar & u);
    void print(const Vector & u, const Comp & m);

    void save(const Scalar & phi, const char * nm, const int it);
    void save(const Vector & u, const Comp & m,
              const char * nm, const int it);

    void save_grid(const Scalar & phi, const char * nm);
    void save_grid(const Vector & u, const Comp & m,
              const char * nm);

  private:
    const char   * name;
    std::vector<Location *> mons; 
    Range<int> r_i;
    Range<int> r_j;
    Range<int> r_k;
};

#endif
