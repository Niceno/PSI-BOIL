/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to/in Grid1D).
*
*  Used for setting (and checking) the location of bnd buffer cells in Grid1D.
*******************************************************************************/
#ifndef BNDGRID_H
#define BNDGRID_H

///////////////
//           //
//  BndGrid  //
//           //
///////////////
class BndGrid {
  public:
    static const BndGrid undefined()   {return BndGrid(0);}
    static const BndGrid wall()        {return BndGrid(1);}
    static const BndGrid extrapolate() {return BndGrid(2);}
    static const BndGrid symmetry()    {return BndGrid(3);}
    
    operator int () const {return val;}
    
  private:
    /* prevent creation of new cutoff conditions */
    explicit BndGrid(const int i) {val = i;}
    int val;
};

#endif
