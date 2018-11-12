/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to/in Grid1D).
*
*  Used for setting (and checking) the periodicity in Grid1D.
*******************************************************************************/
#ifndef PERIODIC_H
#define PERIODIC_H

////////////////
//            //
//  Periodic  //
//            //
////////////////
class Periodic {
  public:
    static const Periodic yes() {return Periodic(true);}
    static const Periodic no()  {return Periodic(false);}
    
    operator bool () const {return val;}
    
  private:
    /* prevent creation of new periodicity conditions */
    explicit Periodic(const bool i) {val = i;}
    bool val;
};

#endif
