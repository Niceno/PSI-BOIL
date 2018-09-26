#ifndef COMP_VECTOR_H
#define COMP_VECTOR_H

#include <iostream>
#include <cassert>

////////////
//        //
//  Comp  //
//        //
////////////
class Comp {

  public:
    Comp()            {val=-1;}

    static const Comp undefined()   {return Comp(-1);}
    static const Comp i()           {return Comp( 0);}
    static const Comp j()           {return Comp( 1);}
    static const Comp k()           {return Comp( 2);}
    static const Comp u()           {return i();}
    static const Comp v()           {return j();}
    static const Comp w()           {return k();}
    static const Comp inf()         {return Comp( 0);}
    static const Comp coefficient() {return Comp( 1);}

    //! Prints the components name.
    friend std::ostream &
      operator << (std::ostream & os, const Comp & com);

    //! Operators.
    Comp & operator ++ ()  // prefix
     {val++; return *this;}
    Comp operator ++ (int) // postfix
     {Comp before(*this); val++; return before;}
    Comp & operator -- ()  // prefix
     {val--; return *this;}
    Comp operator -- (int) // postfix
     {Comp before(*this); val--; return before;}
    bool operator == (const Comp & o) const {return val == o.val;}
    bool operator != (const Comp & o) const {return val != o.val;}
    bool operator <  (const Comp & o) const {return val <  o.val;}
    bool operator >  (const Comp & o) const {return val >  o.val;}
    bool operator <= (const Comp & o) const {return val <= o.val;}
    bool operator >= (const Comp & o) const {return val >= o.val;}

    int operator ~ () const {return val;}
//  operator int () {return val;}

  private:
    int val;

    /* avoid implicit conversions of integer to Comp */
    explicit Comp(const int m) {val = m;} 
};

#endif

/*-----------------------------------------------------------------------------+
 '$Id: comp.h,v 1.3 2016/02/26 11:07:43 niceno Exp $'/
+-----------------------------------------------------------------------------*/
