#ifndef ACCURACY_ORDER_H
#define ACCURACY_ORDER_H

#include <cstdlib>
#include <iostream>
#include <fstream>

/////////////////////
//                 //
// Accuracy order  //
//                 //
/////////////////////
/* this is a ravioli class for accuracy order selection */
class AccuracyOrder {
  public:
    AccuracyOrder() {val=-1;}
    explicit AccuracyOrder(const int m, const bool upwind = false) {
      val = m;
      upw = upwind; /* for the moment, only affects first order */
    }

    static const AccuracyOrder undefined() {return AccuracyOrder(-1);}
    static const AccuracyOrder Zeroth()    {return AccuracyOrder( 0);}
    static const AccuracyOrder First()     {return AccuracyOrder( 1);}
    static const AccuracyOrder FirstCentral()  {return AccuracyOrder(1,false);}
    static const AccuracyOrder FirstUpwind()   {return AccuracyOrder(1, true);}
    static const AccuracyOrder Second()    {return AccuracyOrder( 2);}
    static const AccuracyOrder Third()     {return AccuracyOrder( 3);}
    static const AccuracyOrder Fourth()    {return AccuracyOrder( 4);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const AccuracyOrder & com) {
      ost << com.val<<" "<<com.upw;
      return ost;
    }

    bool operator == (const AccuracyOrder & o) const 
      { return (val == o.val) && (upw == o.upw); }
    bool operator != (const AccuracyOrder & o) const 
      { return (val != o.val) && (upw != o.upw); }

    int eval() const { return val; }
    bool upwind() const { return upw; }

  private:
    int val;
    bool upw;

};
#endif
