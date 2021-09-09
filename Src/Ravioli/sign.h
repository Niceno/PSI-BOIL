#ifndef SIGN_H
#define SIGN_H

#include <cstdlib>
#include <iostream>
#include <fstream>

/***************************************************************************//**
*  \brief Ravioli class for safe specification of positive/negative direction.
*
*  Represents the direction +/-.
*******************************************************************************/

////////////
//        //
//  Sign  //
//        //
////////////
class Sign {
  public:
    Sign() {val=0;} /* undefined */

    static const Sign undefined()   {return Sign(0);}
    static const Sign pos() {return Sign(+1);}
    static const Sign neg() {return Sign(-1);}
    
    //! Prints the sign name.
    friend std::ostream & 
      operator << (std::ostream & os, const Sign & sig);

    //! Operators.
    Sign operator - () const {return Sign(-val);}
    operator int () const {return val;}
    Sign & operator *= (const Sign & s) {
      val *= int(s);
      return *this;
    }
    
    friend Sign operator * (Sign first, const Sign & second) {
      first *= second;
      return first;
    }

  private:
    /* prevent creation of new signs */
    explicit Sign(const int i) {val = i;}
    int val;
};
#endif
