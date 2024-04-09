#ifndef RANGE_H
#define RANGE_H

#include <iostream>

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
    bool    contains(const T i) const {return ( (i>=a) && (i<=b) );}
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
