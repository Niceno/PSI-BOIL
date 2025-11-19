#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include "../Global/global_precision.h"

template <class T> class Range;
template <class T> std::ostream & operator << (std::ostream &, 
                                               const Range<T> &);
/////////////
//         //
//  Range  //
//         //
/////////////
template <class T> class Range {
  public:
    explicit Range(const T i, const T j) : a(i), b(j) {}
    explicit Range() {a=0; b=-1;}
    const T first() const {return a;}
    const T last () const {return b;}
    void    first(const T i) {a=i;}
    void    last (const T j) {b=j;}
    const T dist () const { return b-a; }
    bool    contains(const T i) const {return ( (i>=a) && (i<=b) );}
    real fraction(const T i) const {
      return std::max(0.0,
               std::min(1.0,
                 (i-a)/(b-a)));
    }
    bool    exists() const {return a<=b;}

    friend std::ostream & 
      operator << <> (std::ostream & os, const Range<T> & rng);

  private:
    T a, b;
};

/******************************************************************************/
template <class T> std::ostream & operator<< (std::ostream &ost,
                                              const Range<T> & rng) {
  ost << rng.first() << " -> " << rng.last ();

  return ost;
}

#endif
