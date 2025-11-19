#ifndef DIR_H
#define DIR_H

#include <cstdlib>
#include <iostream>
#include <fstream>

/***************************************************************************//**
*  \brief Ravioli class for safe specification of logical (i,j,k) direction.
*
*  Represents the direction in logical coordinates (i,j,k). Used extensivelly
*  in many parts of the code for setting boundary conditions, periodic 
*  boundaries, etc.
*******************************************************************************/

///////////
//       //
//  Dir  //
//       //
///////////
class Dir {
  public:
    static const Dir undefined() {return Dir(-1);}
    static const Dir imin()      {return Dir( 0);}
    static const Dir imax()      {return Dir( 1);}
    static const Dir jmin()      {return Dir( 2);}
    static const Dir jmax()      {return Dir( 3);}
    static const Dir kmin()      {return Dir( 4);}
    static const Dir kmax()      {return Dir( 5);}
    static const Dir ibody()     {return Dir( 6);}
    
    operator int () const {return val;}
    
    //! Prints the direction name.
    friend std::ostream & 
      operator << (std::ostream & os, const Dir & dir);

  private:
    /* prevent creation of new b.c. */
    explicit Dir(const int i) {val = i;}
    int val;
};

#endif
