/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to Plot).
*
*  Used for setting whether the buffers will be plotted.
*******************************************************************************/
#ifndef BUFFERS_H
#define BUFFERS_H

///////////////
//           //
//  Buffers  //
//           //
///////////////
class Buffers {
  public:
    static const Buffers yes() {return Buffers(true);}
    static const Buffers no () {return Buffers(false);}
    
    operator bool () const {return val;}
    
  private:
    /* prevent creation of new decompose conditions */
    explicit Buffers(const bool i) {val = i;}
    bool val;
};	

#endif
