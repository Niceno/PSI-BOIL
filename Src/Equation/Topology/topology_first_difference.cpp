#include "topology.h"

/* 
   presented expressions were calculated using symbolic python package sympy,
   based on theory from 
     https://en.wikipedia.org/wiki/Finite_difference_coefficient#cite_note-4
     http://web.media.mit.edu/~crtaylor/calculator.html
*/

/******************************************************************************/
real Topology::nth_order_first(const std::vector<StencilPoint> & stencil,
                               const AccuracyOrder & order) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a desired-order difference.
*******************************************************************************/
  switch(order.eval()) {
    case 0 :
      return zeroth_order_first(stencil);
    case 1 :
      return first_order_first(stencil);
    case 2 :
      return second_order_first(stencil);
    case 3 :
      return third_order_first(stencil);
    case 4 :
      return fourth_order_first(stencil);
    default :
      boil::aout<<"Topology: unrecognised first difference requested. Exiting."
                <<boil::endl;
      exit(0);
  }

  return 0.0;
}

/******************************************************************************/
real Topology::zeroth_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a zeroth-order difference.
*******************************************************************************/
  return 0.0;
}

/******************************************************************************/
real Topology::first_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a first-order difference.
*******************************************************************************/
  std::vector<real> coefs;
  first_order_first_coefs(coefs,stencil);

  return coefs[0]*stencil[0].val+coefs[1]*stencil[1].val;

}

/******************************************************************************/
real Topology::second_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  std::vector<real> coefs;
  second_order_first_coefs(coefs,stencil);

  return coefs[0]*stencil[0].val+coefs[1]*stencil[1].val
        +coefs[2]*stencil[2].val;               
}

/******************************************************************************/
real Topology::third_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a third-order difference.
*******************************************************************************/
  std::vector<real> coefs;
  third_order_first_coefs(coefs,stencil);

  return coefs[0]*stencil[0].val+coefs[1]*stencil[1].val
        +coefs[2]*stencil[2].val+coefs[3]*stencil[3].val;
}

/******************************************************************************/
real Topology::fourth_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a fourth-order difference.
*******************************************************************************/
  std::vector<real> coefs;
  fourth_order_first_coefs(coefs,stencil);
  
  return coefs[0]*stencil[0].val+coefs[1]*stencil[1].val
        +coefs[2]*stencil[2].val+coefs[3]*stencil[3].val
        +coefs[4]*stencil[4].val;
}
