/***************************************************************************//**
*  \brief Ravioli class for safer parameter passing (to Plot).
*
*  Used for setting whether the cell-centered data will be plotted as nodal.
*******************************************************************************/
#ifndef ASNODES_H
#define ASNODES_H

///////////////
//           //
//  AsNodes  //
//           //
///////////////
class AsNodes {
  public:
    static const AsNodes yes() {return AsNodes(true);}
    static const AsNodes no()  {return AsNodes(false);}
    
    operator bool () const {return val;}
    
  private:
    /* prevent creation of new decompose conditions */
    explicit AsNodes(const bool i) {val = i;}
    bool val;
};	

#endif

/*-----------------------------------------------------------------------------+
 '$Id: asnodes.h,v 1.1 2009/10/29 14:20:58 niceno Exp $'/
+-----------------------------------------------------------------------------*/
