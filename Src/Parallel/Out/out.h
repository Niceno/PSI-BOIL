#ifndef PAROUT_H
#define PAROUT_H

#include "../communicator.h"

#include <iostream>
#include <sstream>

/***************************************************************************//**
*  \brief Ravioli class used for parallel output streams.       
*
*  Empty class, only instructing parallel output stream AllOut to print 
*  processor i.d.  
*******************************************************************************/

//////////////
//          //
//  ProcID  //
//          //
//////////////
class ProcID {};

/***************************************************************************//**
*  \brief Parallel output stream. All processors print.
*
*  Used in occasions when each processor has to print data to screen. That is
*  usually (but not always) needed for debugging. Global object "aout", 
*  defined in namespace "boil" is created to ease the use of it. It is used
*  as:
*
*  \code
*    boil::aout << "Hi from a processor!" << boil::endl;
*  \endcode
*******************************************************************************/

//////////////
//          //
//  AllOut  //
//          //
//////////////
class AllOut : public std::ostream {
  public:
    AllOut() : std::ostream(std::cout.rdbuf()) {}

    AllOut & operator << (const char * c) 
     {std::cout << c; return *this;} 

    AllOut & operator << (const char c) 
     {std::cout << c; return *this;} 

    AllOut & operator << (const int i) 
     {std::cout << i; return *this;} 

    AllOut & operator << (const unsigned u) 
     {std::cout << u; return *this;} 

    AllOut & operator << (const real & r) 
     {std::cout << r; return *this;} 

    AllOut & operator << (const std::string & s) 
     {std::cout << s; return *this;} 

    AllOut & operator << (const std::stringstream & s) 
     {std::cout << s.str(); return *this;} 

    AllOut & operator << (const ProcID & p) 
     {std::cout << boil::cart.iam() << ": "; return *this;} 

    AllOut & operator << ( std::ostream & (*f) (std::ostream&) ) 
     {f(*this); return *this;} 

    AllOut & operator << (const Dir & d) 
     {std::cout << d; return *this;} 
};

/***************************************************************************//**
*  \brief Parallel output stream. Only one processor prints.
*
*  Used in occasions when only one processor has to print data to screen. 
*  That is what is usually needed in a parallel application. Global object 
*  "oout", defined in namespace "boil", is created to ease the use of it. 
*  It is used as:
*
*  \code
*    boil::oout << "Start of simulation" << boil::endl;
*  \endcode
*******************************************************************************/

//////////////
//          //
//  OneOut  //
//          //
//////////////
class OneOut : public std::ostream {
  public:
    OneOut() : std::ostream(std::cout.rdbuf()) {}

    OneOut & operator << (const char * c) 
     {if(!boil::cart.iam()) std::cout << c; return *this;} 

    OneOut & operator << (const char c) 
     {if(!boil::cart.iam()) std::cout << c; return *this;} 

    OneOut & operator << (const int i) 
     {if(!boil::cart.iam()) std::cout << i; return *this;} 

    OneOut & operator << (const unsigned u) 
     {std::cout << u; return *this;} 

    OneOut & operator << (const real & r) 
     {if(!boil::cart.iam()) std::cout << r; return *this;} 

    OneOut & operator << (const std::string & s) 
     {if(!boil::cart.iam()) std::cout << s; return *this;} 

    OneOut & operator << (const std::stringstream & s) 
     {if(!boil::cart.iam()) std::cout << s.str(); return *this;} 

    OneOut & operator << ( std::ostream & (*f) (std::ostream&) ) 
     {if(!boil::cart.iam()) f(*this); return *this;} 

    OneOut & operator << (const Dir & d) 
     {if(!boil::cart.iam()) std::cout << d; return *this;} 
};

/*-----------------+
|  global streams  |
+-----------------*/
namespace boil {
  extern AllOut aout;
  extern OneOut oout;
  extern ProcID pid;
}

#endif

/*-----------------------------------------------------------------------------+
 '$Id: out.h,v 1.12 2017/10/20 12:23:19 sato Exp $'/
+-----------------------------------------------------------------------------*/
