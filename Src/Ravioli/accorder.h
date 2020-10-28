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
    explicit AccuracyOrder(const int m) {val = m;}

    static const AccuracyOrder undefined() {return AccuracyOrder(-1);}
    static const AccuracyOrder Zeroth()    {return AccuracyOrder( 0);}
    static const AccuracyOrder First()     {return AccuracyOrder( 1);}
    static const AccuracyOrder Second()    {return AccuracyOrder( 2);}
    static const AccuracyOrder Third()     {return AccuracyOrder( 3);}
    static const AccuracyOrder Fourth()    {return AccuracyOrder( 4);}

    //! Prints the components name.
    friend std::ostream & operator << (std::ostream & ost, const AccuracyOrder & com) {
      ost << com.val;
      return ost;
    }

    bool operator == (const AccuracyOrder & o) const {return val == o.val;}
    bool operator != (const AccuracyOrder & o) const {return val != o.val;}

    int eval() const { return val; }

  private:
    int val;

};
#endif
