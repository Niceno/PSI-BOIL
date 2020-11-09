#include "topology.h"

/* 
   presented expressions were calculated using symbolic python package sympy,
   based on theory from 
     https://en.wikipedia.org/wiki/Finite_difference_coefficient#cite_note-4
     http://web.media.mit.edu/~crtaylor/calculator.html
*/

/******************************************************************************/
real Topology::nth_order_third(const std::vector<StencilPoint> & stencil,
                               const AccuracyOrder & order) const {
/***************************************************************************//**
*  \brief Approximate third derivative using a desired-order difference.
*******************************************************************************/
  switch(order.eval()) {
    case 0 :
      return zeroth_order_third(stencil);
    case 1 :
      return first_order_third(stencil);
    case 2 :
      return third_order_third(stencil);
    case 3 :
      return third_order_third(stencil);
    case 4 :
      return fourth_order_third(stencil);
    default :
      boil::aout<<"Topology: unrecognised third difference requested. Exiting."
                <<boil::endl;
      exit(0);
  }

  return 0.0;
}

/******************************************************************************/
real Topology::zeroth_order_third(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate third derivative using a zeroth-order difference.
*******************************************************************************/
  return 0.0;
}

/******************************************************************************/
real Topology::first_order_third(const std::vector<StencilPoint> & stencil)
                                                                       const {
/***************************************************************************//**
*  \brief Approximate third derivative using a first-order difference.
*******************************************************************************/
  return 0.0;
}

/******************************************************************************/
real Topology::second_order_third(const std::vector<StencilPoint> & stencil)
                                                                       const {
/***************************************************************************//**
*  \brief Approximate third derivative using a second-order difference.
*******************************************************************************/
  return 0.0;
}

/******************************************************************************/
real Topology::third_order_third(const std::vector<StencilPoint> & stencil)
                                                                       const {
/***************************************************************************//**
*  \brief Approximate third derivative using a third-order difference.
*******************************************************************************/

  real c0 =  6./(stencil[0].pos*stencil[0].pos*stencil[0].pos
               - stencil[0].pos*stencil[0].pos*stencil[1].pos
               - stencil[0].pos*stencil[0].pos*stencil[2].pos
               - stencil[0].pos*stencil[0].pos*stencil[3].pos
               + stencil[0].pos*stencil[1].pos*stencil[2].pos
               + stencil[0].pos*stencil[1].pos*stencil[3].pos
               + stencil[0].pos*stencil[2].pos*stencil[3].pos
               - stencil[1].pos*stencil[2].pos*stencil[3].pos);

  real c1 = -6./(stencil[0].pos*stencil[1].pos*stencil[1].pos
               - stencil[0].pos*stencil[1].pos*stencil[2].pos
               - stencil[0].pos*stencil[1].pos*stencil[3].pos
               + stencil[0].pos*stencil[2].pos*stencil[3].pos
               - stencil[1].pos*stencil[1].pos*stencil[1].pos
               + stencil[1].pos*stencil[1].pos*stencil[2].pos
               + stencil[1].pos*stencil[1].pos*stencil[3].pos
               - stencil[1].pos*stencil[2].pos*stencil[3].pos);
  
  real c2 =  6./(stencil[0].pos*stencil[1].pos*stencil[2].pos
               - stencil[0].pos*stencil[1].pos*stencil[3].pos
               - stencil[0].pos*stencil[2].pos*stencil[2].pos
               + stencil[0].pos*stencil[2].pos*stencil[3].pos
               - stencil[1].pos*stencil[2].pos*stencil[2].pos
               + stencil[1].pos*stencil[2].pos*stencil[3].pos
               + stencil[2].pos*stencil[2].pos*stencil[2].pos
               - stencil[2].pos*stencil[2].pos*stencil[3].pos);

  real c3 = -6./(stencil[0].pos*stencil[1].pos*stencil[2].pos
               - stencil[0].pos*stencil[1].pos*stencil[3].pos
               - stencil[0].pos*stencil[2].pos*stencil[3].pos
               + stencil[0].pos*stencil[3].pos*stencil[3].pos 
               - stencil[1].pos*stencil[2].pos*stencil[3].pos
               + stencil[1].pos*stencil[3].pos*stencil[3].pos
               + stencil[2].pos*stencil[3].pos*stencil[3].pos 
               - stencil[3].pos*stencil[3].pos*stencil[3].pos);

  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val;
}

/******************************************************************************/
real Topology::fourth_order_third(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate third derivative using a fourth-order difference.
*******************************************************************************/

  real c0 =  -6.*(stencil[1].pos + stencil[2].pos
                + stencil[3].pos + stencil[4].pos)
                /(stencil[0].pos*stencil[0].pos*stencil[0].pos*stencil[0].pos
                - stencil[0].pos*stencil[0].pos*stencil[0].pos*stencil[1].pos
                - stencil[0].pos*stencil[0].pos*stencil[0].pos*stencil[2].pos
                - stencil[0].pos*stencil[0].pos*stencil[0].pos*stencil[3].pos
                - stencil[0].pos*stencil[0].pos*stencil[0].pos*stencil[4].pos
                + stencil[0].pos*stencil[0].pos*stencil[1].pos*stencil[2].pos
                + stencil[0].pos*stencil[0].pos*stencil[1].pos*stencil[3].pos
                + stencil[0].pos*stencil[0].pos*stencil[1].pos*stencil[4].pos
                + stencil[0].pos*stencil[0].pos*stencil[2].pos*stencil[3].pos
                + stencil[0].pos*stencil[0].pos*stencil[2].pos*stencil[4].pos
                + stencil[0].pos*stencil[0].pos*stencil[3].pos*stencil[4].pos
                - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
                - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
                - stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
                - stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos 
                + stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos);

  real c1 =  6.*(stencil[0].pos + stencil[2].pos
               + stencil[3].pos + stencil[4].pos)
               /(stencil[0].pos*stencil[1].pos*stencil[1].pos*stencil[1].pos
               - stencil[0].pos*stencil[1].pos*stencil[1].pos*stencil[2].pos
               - stencil[0].pos*stencil[1].pos*stencil[1].pos*stencil[3].pos
               - stencil[0].pos*stencil[1].pos*stencil[1].pos*stencil[4].pos
               + stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
               + stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
               + stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
               - stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               - stencil[1].pos*stencil[1].pos*stencil[1].pos*stencil[1].pos
               + stencil[1].pos*stencil[1].pos*stencil[1].pos*stencil[2].pos
               + stencil[1].pos*stencil[1].pos*stencil[1].pos*stencil[3].pos
               + stencil[1].pos*stencil[1].pos*stencil[1].pos*stencil[4].pos
               - stencil[1].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
               - stencil[1].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
               - stencil[1].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
               + stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos);

  real c2 = -6.*(stencil[0].pos + stencil[1].pos
               + stencil[3].pos + stencil[4].pos)
               /(stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[2].pos
               - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
               - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
               + stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
               - stencil[0].pos*stencil[2].pos*stencil[2].pos*stencil[2].pos
               + stencil[0].pos*stencil[2].pos*stencil[2].pos*stencil[3].pos
               + stencil[0].pos*stencil[2].pos*stencil[2].pos*stencil[4].pos
               - stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               - stencil[1].pos*stencil[2].pos*stencil[2].pos*stencil[2].pos
               + stencil[1].pos*stencil[2].pos*stencil[2].pos*stencil[3].pos
               + stencil[1].pos*stencil[2].pos*stencil[2].pos*stencil[4].pos
               - stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               + stencil[2].pos*stencil[2].pos*stencil[2].pos*stencil[2].pos
               - stencil[2].pos*stencil[2].pos*stencil[2].pos*stencil[3].pos
               - stencil[2].pos*stencil[2].pos*stencil[2].pos*stencil[4].pos
               + stencil[2].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos);

  real c3 =  6.*(stencil[0].pos + stencil[1].pos
               + stencil[2].pos + stencil[4].pos)
               /(stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
               - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
               - stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[3].pos
               + stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
               - stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[3].pos
               + stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               + stencil[0].pos*stencil[3].pos*stencil[3].pos*stencil[3].pos
               - stencil[0].pos*stencil[3].pos*stencil[3].pos*stencil[4].pos
               - stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[3].pos
               + stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               + stencil[1].pos*stencil[3].pos*stencil[3].pos*stencil[3].pos
               - stencil[1].pos*stencil[3].pos*stencil[3].pos*stencil[4].pos
               + stencil[2].pos*stencil[3].pos*stencil[3].pos*stencil[3].pos
               - stencil[2].pos*stencil[3].pos*stencil[3].pos*stencil[4].pos
               - stencil[3].pos*stencil[3].pos*stencil[3].pos*stencil[3].pos
               + stencil[3].pos*stencil[3].pos*stencil[3].pos*stencil[4].pos);

  real c4 = -6.*(stencil[0].pos + stencil[1].pos
               + stencil[2].pos + stencil[3].pos)
               /(stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[3].pos
               - stencil[0].pos*stencil[1].pos*stencil[2].pos*stencil[4].pos
               - stencil[0].pos*stencil[1].pos*stencil[3].pos*stencil[4].pos
               + stencil[0].pos*stencil[1].pos*stencil[4].pos*stencil[4].pos
               - stencil[0].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               + stencil[0].pos*stencil[2].pos*stencil[4].pos*stencil[4].pos
               + stencil[0].pos*stencil[3].pos*stencil[4].pos*stencil[4].pos
               - stencil[0].pos*stencil[4].pos*stencil[4].pos*stencil[4].pos
               - stencil[1].pos*stencil[2].pos*stencil[3].pos*stencil[4].pos
               + stencil[1].pos*stencil[2].pos*stencil[4].pos*stencil[4].pos
               + stencil[1].pos*stencil[3].pos*stencil[4].pos*stencil[4].pos
               - stencil[1].pos*stencil[4].pos*stencil[4].pos*stencil[4].pos
               + stencil[2].pos*stencil[3].pos*stencil[4].pos*stencil[4].pos
               - stencil[2].pos*stencil[4].pos*stencil[4].pos*stencil[4].pos
               - stencil[3].pos*stencil[4].pos*stencil[4].pos*stencil[4].pos
               + stencil[4].pos*stencil[4].pos*stencil[4].pos*stencil[4].pos);

  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val+c4*stencil[4].val;
}
