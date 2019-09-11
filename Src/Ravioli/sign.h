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
    Sign()            {val=-1;}

    static const Sign undefined()   {return Sign(0);}
    static const Sign pos() {return Sign(+1);}
    static const Sign neg() {return Sign(-1);}
    
    operator int () const {return val;}
    
    //! Prints the sign name.
    friend std::ostream & 
      operator << (std::ostream & os, const Sign & sig);

  private:
    /* prevent creation of new signs */
    explicit Sign(const int i) {val = i;}
    int val;
};

#endif
