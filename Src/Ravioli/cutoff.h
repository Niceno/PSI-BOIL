/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to/in Grid1D).
*
*  Used for setting (and checking) the location of bnd buffer cells in Grid1D.
*******************************************************************************/
#ifndef CUTOFF_H
#define CUTOFF_H

//////////////
//          //
//  Cutoff  //
//          //
//////////////
class Cutoff {
  public:
    static const Cutoff undefined()   {return Cutoff(0);}
    static const Cutoff wall()        {return Cutoff(1);}
    static const Cutoff extrapolate() {return Cutoff(2);}
    static const Cutoff symmetry()    {return Cutoff(3);}
    
    operator int () const {return val;}
    
  private:
    /* prevent creation of new cutoff conditions */
    explicit Cutoff(const int i) {val = i;}
    int val;
};

#endif
