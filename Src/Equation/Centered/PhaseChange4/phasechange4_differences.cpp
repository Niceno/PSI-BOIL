#include "phasechange4.h"

/* 
   the presented expression were calculated using symbolic python package sympy,
   based on theory from 
     https://en.wikipedia.org/wiki/Finite_difference_coefficient#cite_note-4
     http://web.media.mit.edu/~crtaylor/calculator.html
*/

/******************************************************************************/
real PhaseChange4::second_order_difference(const std::vector<real> & stencil,
                                           const std::vector<real> & values) {
/***************************************************************************//**
*  \brief Approximate first order derivative using a second-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1] - 1./stencil[2];
  real c1 = -stencil[2]/(stencil[1]*(stencil[1]-stencil[2]));
  real c2 =  stencil[1]/(stencil[2]*(stencil[1]-stencil[2]));

  return c0*values[0]+c1*values[1]+c2*values[2];
}

/******************************************************************************/
real PhaseChange4::third_order_difference(const std::vector<real> & stencil,
                                          const std::vector<real> & values) {
/***************************************************************************//**
*  \brief Approximate first order derivative using a third-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1] - 1./stencil[2] - 1./stencil[3]; 
  
  real c1 =  stencil[2]*stencil[3]
           /(stencil[1]*(stencil[1]*stencil[1]
                       - stencil[1]*stencil[2]
                       - stencil[1]*stencil[3]
                       + stencil[2]*stencil[3]));

  real c2 = -stencil[1]*stencil[3]
           /(stencil[2]*(stencil[1]*stencil[2]
                       - stencil[1]*stencil[3]
                       - stencil[2]*stencil[2]
                       + stencil[2]*stencil[3]));

  real c3 =  stencil[1]*stencil[2]
           /(stencil[3]*(stencil[1]*stencil[2]
                       - stencil[1]*stencil[3]
                       - stencil[2]*stencil[3]
                       + stencil[3]*stencil[3]));

  return c0*values[0]+c1*values[1]+c2*values[2]+c3*values[3];
}

/******************************************************************************/
real PhaseChange4::fourth_order_difference(const std::vector<real> & stencil,
                                           const std::vector<real> & values) {
/***************************************************************************//**
*  \brief Approximate first order derivative using a fourth-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1] - 1./stencil[2] - 1./stencil[3] - 1./stencil[4]; 

  real c1 = -stencil[2]*stencil[3]*stencil[4]
           /(stencil[1]*(stencil[1]*stencil[1]*stencil[1]
                       - stencil[1]*stencil[1]*stencil[2]
                       - stencil[1]*stencil[1]*stencil[3]
                       - stencil[1]*stencil[1]*stencil[4]
                       + stencil[1]*stencil[2]*stencil[3]
                       + stencil[1]*stencil[2]*stencil[4]
                       + stencil[1]*stencil[3]*stencil[4]
                       - stencil[2]*stencil[3]*stencil[4]));

  real c2 =  stencil[1]*stencil[3]*stencil[4]
           /(stencil[2]*(stencil[1]*stencil[2]*stencil[2]
                       - stencil[1]*stencil[2]*stencil[3]
                       - stencil[1]*stencil[2]*stencil[4]
                       + stencil[1]*stencil[3]*stencil[4]
                       - stencil[2]*stencil[2]*stencil[2]
                       + stencil[2]*stencil[2]*stencil[3]
                       + stencil[2]*stencil[2]*stencil[4]
                       - stencil[2]*stencil[3]*stencil[4]));

  real c3 = -stencil[1]*stencil[2]*stencil[4]
           /(stencil[3]*(stencil[1]*stencil[2]*stencil[3]
                       - stencil[1]*stencil[2]*stencil[4]
                       - stencil[1]*stencil[3]*stencil[3]
                       + stencil[1]*stencil[3]*stencil[4]
                       - stencil[2]*stencil[3]*stencil[3]
                       + stencil[2]*stencil[3]*stencil[4]
                       + stencil[3]*stencil[3]*stencil[3]
                       - stencil[3]*stencil[3]*stencil[4]));

  real c4 =  stencil[1]*stencil[2]*stencil[3]
           /(stencil[4]*(stencil[1]*stencil[2]*stencil[3]
                       - stencil[1]*stencil[2]*stencil[4]
                       - stencil[1]*stencil[3]*stencil[4]
                       + stencil[1]*stencil[4]*stencil[4]
                       - stencil[2]*stencil[3]*stencil[4]
                       + stencil[2]*stencil[4]*stencil[4]
                       + stencil[3]*stencil[4]*stencil[4]
                       - stencil[4]*stencil[4]*stencil[4]));
  
  return c0*values[0]+c1*values[1]+c2*values[2]+c3*values[3]+c4*values[4];
}
