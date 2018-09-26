#include "global_limiter.h"

/***************************************************************************//**
*  Prints the name of the current limiter. Needed only for information or 
*  debugging.
*******************************************************************************/
std::ostream & operator << (std::ostream &ost, const ConvScheme & cs) {

  switch(cs.val) {
    case( 1): ost << "upwind ";   break;
    case( 2): ost << "central ";  break;
    case( 3): ost << "minmod ";   break;
    case( 4): ost << "smart ";    break;
    case( 5): ost << "muscl ";    break;
    case( 6): ost << "superbee "; break;
    case( 7): ost << "van_leer";  break;
    case( 8): ost << "mc";        break;
  }
  
  return ost;
}

/***************************************************************************//**
*  Sets the convection scheme at the time of creation.                  
*******************************************************************************/
Limiter::Limiter(const ConvScheme & cs) {

       if(cs == ConvScheme::upwind()  ) b = & Limiter::upwind;
  else if(cs == ConvScheme::central() ) b = & Limiter::central;
  else if(cs == ConvScheme::minmod()  ) b = & Limiter::minmod;
  else if(cs == ConvScheme::smart ()  ) b = & Limiter::smart;
  else if(cs == ConvScheme::muscl ()  ) b = & Limiter::muscl;
  else if(cs == ConvScheme::superbee()) b = & Limiter::superbee;
  else if(cs == ConvScheme::van_leer()) b = & Limiter::van_leer;
  else if(cs == ConvScheme::mc()      ) b = & Limiter::mc;
  else {
    boil::oout << "unknown limiter! exiting!" << boil::endl;
    exit(0);
  }
}

/***************************************************************************//**
*  Sets the convection scheme. Uses argument of ravioli type ConvScheme
*  for safe argument passing.
*******************************************************************************/
void Limiter::set(const ConvScheme & cs) {

       if(cs == ConvScheme::upwind()  ) b = & Limiter::upwind;
  else if(cs == ConvScheme::central() ) b = & Limiter::central;
  else if(cs == ConvScheme::minmod()  ) b = & Limiter::minmod;
  else if(cs == ConvScheme::smart()   ) b = & Limiter::smart;
  else if(cs == ConvScheme::muscl()   ) b = & Limiter::muscl;
  else if(cs == ConvScheme::superbee()) b = & Limiter::superbee;
  else if(cs == ConvScheme::van_leer()) b = & Limiter::van_leer;
  else if(cs == ConvScheme::mc()      ) b = & Limiter::mc;
  else {
    boil::oout << "unknown limiter! exiting!" << boil::endl;
    exit(0);
  }
}

/***************************************************************************//**
*
*  \param     ucp  - velocity at the downwind face,
*  \param     phic - value at the cell centre,
*  \param     phim - value at the upwind cell centre,
*  \param     phip - value at the downwind cell centre.
*
*  Implements a generalized limiter. It is a framework algorith for all
*  implemented limiters. Which one is used is determined from the local
*  function pointer "b". 
*
*  This implementation Uses information from thee cell centres
*  to limit the value at the cell face. To avoid the usage of extended 
*  stencils (more than one neighbour in each direction), it compute (limits)
*  only the downwind cell face value. 
*
*  For example, if for a cell at i,j,k (denoted by "c", as central), velocity
*  at the east face (denoted here by "cp") flows to the right, we have all 
*  the necessary values (at "c", "m" and "p") to compute a limiting value. 
*
*  \image   html  limiter_a.jpg "Limiter on the east face"
*  \image   latex limiter_a.eps "Limiter on the easr face" width=10cm
*
*  If the velocity at the east face was headed to the left, the computation
*  of the value on this cell-face would be left for the next cell (i+1,j,k).
*  In the same spirit, if the velocity on the west face was oriented to the
*  left, we have all the values needed for the limiter.
*
*  \image   html  limiter_b.jpg "Limiter on the west face"
*  \image   latex limiter_b.eps "Limiter on the west face" width=10cm
*
*******************************************************************************/
real Limiter::limit(const real ucp,
                    const real phim,
                    const real phic,
                    const real phip) {

  if(ucp > 0.0) {
    /* avoid division by zero */
    if( fabs(phic-phim) > boil::atto ) {
      real r = (phip - phic) / (phic - phim);
      real psi = (*this.*b)(r); 
      return phic + 0.5 * psi * (phic - phim);
    /* / 0 is represented as * 1e+18; * 0 is represented as * 1e-18 */
    } else {
      real r = (phip - phic) * boil::exa; 
      real psi = (*this.*b)(r);
      return phic + 0.5 * psi * boil::atto;            
    }
  }
  else
    return 0.0;
}

/******************************************************************************/
real Limiter::upwind  (const real r) const 
 {return 0.0;}

/******************************************************************************/
real Limiter::central (const real r) const 
 {return r;}

/******************************************************************************/
real Limiter::minmod  (const real r) const 
 {return boil::maxr( 0.0, boil::minr(1.0, r) );}

/******************************************************************************/
real Limiter::smart   (const real r) const 
 {return boil::maxr( 0.0, boil::minr(4.0, 0.75*r+0.25, 2.0*r) );}

/******************************************************************************/
real Limiter::muscl   (const real r) const 
 {return boil::maxr( 0.0, boil::minr(2.0, 0.5*r+0.5, 2.0*r) );}

/******************************************************************************/
real Limiter::superbee(const real r) const 
 {return boil::maxr( 0.0, boil::minr(2.0*r, 1.0), boil::minr(r,2.0) );}

/******************************************************************************/
real Limiter::van_leer(const real r) const 
 {return (r+fabs(r))/(1.0+fabs(r));}

/******************************************************************************/
real Limiter::mc(const real r) const 
 {return boil::maxr( 0.0, boil::minr(2.0*r, 0.5*(1.0+r), 2.0) );}

/*-----------------------------------------------------------------------------+
 '$Id: global_limiter.cpp,v 1.12 2011/05/25 12:36:04 niceno Exp $'/
+-----------------------------------------------------------------------------*/
