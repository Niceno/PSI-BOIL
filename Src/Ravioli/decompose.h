/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to Domain).
*
*  Used for setting whether the domain will be decomposed for MPI version.
*******************************************************************************/
#ifndef DECOMPOSE_H
#define DECOMPOSE_H

/////////////////
//             //
//  Decompose  //
//             //
/////////////////
class Decompose {
  public:
    static const Decompose no()  {return Decompose(0);}
    static const Decompose x()   {return Decompose(1);}
    static const Decompose y()   {return Decompose(2);}
    static const Decompose z()   {return Decompose(3);}
    static const Decompose xy()  {return Decompose(4);}
    static const Decompose yx()  {return Decompose(4);}
    static const Decompose xz()  {return Decompose(5);}
    static const Decompose zx()  {return Decompose(5);}
    static const Decompose yz()  {return Decompose(6);}
    static const Decompose zy()  {return Decompose(6);}
    static const Decompose xyz() {return Decompose(7);}
    
    operator int () const {return val;}
    
  private:
    /* prevent creation of new decompose conditions */
    explicit Decompose(const int i) {val = i;}
    int val;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: decompose.h,v 1.4 2009/10/27 09:09:47 niceno Exp $'/
+-----------------------------------------------------------------------------*/
