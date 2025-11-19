#ifndef BOX_H
#define BOX_H

#include <iostream>

#include "range.h"
#include "comp.h"

template <class T> class Box;
template <class T> std::ostream & operator << (std::ostream &, 
                                               const Box<T> &);
///////////
//       //
//  Box  //
//       //
///////////
template <class T> class Box {
  public:
    explicit Box(const Range<T> & i, const Range<T> & j, const Range<T> & k) 
     {a[0]=i; a[1]=j; a[2]=k;}
    explicit Box() {}
    const T first(const Comp & c) const {return a[~c].first();}
    const T last (const Comp & c) const {return a[~c].last ();}
    void    first(const Comp & c, const T & v) {a[~c].first(v);}
    void    last (const Comp & c, const T & v) {a[~c].last (v);}
    void    move (const Comp & c, const T & v) {a[~c].first(a[~c].first()+v); 
                                                a[~c].last (a[~c].last ()+v);}
    bool contains(const T i, const T j, const T k) const 
     {return a[0].contains(i) && a[1].contains(j) && a[2].contains(k);}
    bool exists() const 
     {return a[0].exists() && a[1].exists() && a[2].exists();}

    friend std::ostream & 
      operator << <> (std::ostream & os, const Box<T> & rng);

  private:
    Range<T> a[3];
};

/******************************************************************************/
template <class T> std::ostream & operator<< (std::ostream &ost,
                                              const Box<T> & rng) {
//  ost << "first = " << rng.first() << ", " 
//      << "last = "  << rng.last ();
  ost << "i: " << rng.a[0].first() << "->" << rng.a[0].last () << ",  "         
      << "j: " << rng.a[1].first() << "->" << rng.a[1].last () << ",  "           
      << "k: " << rng.a[2].first() << "->" << rng.a[2].last ();

  return ost;
}

#endif
