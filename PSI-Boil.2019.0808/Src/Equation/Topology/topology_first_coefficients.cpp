#include "topology.h"

/* 
   presented expressions were calculated using symbolic python package sympy,
   based on theory from 
     https://en.wikipedia.org/wiki/Finite_difference_coefficient#cite_note-4
     http://web.media.mit.edu/~crtaylor/calculator.html
*/

/******************************************************************************/
void Topology::nth_order_first_coefs(std::vector<real> & coefs,
                                     const std::vector<StencilPoint> & stencil,
                                     const AccuracyOrder & order) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a desired-order difference.
*******************************************************************************/
  switch(order.eval()) {
    case 0 :
      return zeroth_order_first_coefs(coefs,stencil);
    case 1 :
      return first_order_first_coefs(coefs,stencil);
    case 2 :
      return second_order_first_coefs(coefs,stencil);
    case 3 :
      return third_order_first_coefs(coefs,stencil);
    case 4 :
      return fourth_order_first_coefs(coefs,stencil);
    default :
      boil::aout<<"Topology: unrecognised first difference requested. Exiting."
                <<boil::endl;
      exit(0);
  }
}

/******************************************************************************/
void Topology::zeroth_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a zeroth-order difference.
*******************************************************************************/
  coefs = { 0. };
  return;
}

/******************************************************************************/
void Topology::first_order_first_coefs(std::vector<real> & coefs,
                                 const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a first-order difference.
*******************************************************************************/
  real c0 =  1./(stencil[0].pos-stencil[1].pos);

  real c1 = -1./(stencil[0].pos-stencil[1].pos);

  coefs = { c0, c1 };
  return;
}

#if 1
/******************************************************************************/
void Topology::second_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  real c0 = -(stencil[1].pos + stencil[2].pos)
            /(stencil[0].pos*stencil[0].pos
            - stencil[0].pos*stencil[1].pos
            - stencil[0].pos*stencil[2].pos
            + stencil[1].pos*stencil[2].pos);

  real c1 =  (stencil[0].pos + stencil[2].pos)
            /(stencil[0].pos*stencil[1].pos
            - stencil[0].pos*stencil[2].pos
            - stencil[1].pos*stencil[1].pos
            + stencil[1].pos*stencil[2].pos);

  real c2 = -(stencil[0].pos + stencil[1].pos)
            /(stencil[0].pos*stencil[1].pos
            - stencil[0].pos*stencil[2].pos
            - stencil[1].pos*stencil[2].pos
            + stencil[2].pos*stencil[2].pos);

  coefs = { c0, c1, c2 };
  return;
}

/******************************************************************************/
void Topology::third_order_first_coefs(std::vector<real> & coefs,
                                 const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a third-order difference.
*******************************************************************************/
  real c0 =  (stencil[1].pos*stencil[2].pos 
            + stencil[1].pos*stencil[3].pos
            + stencil[2].pos*stencil[3].pos)
            /(stencil[0].pos*stencil[0].pos*stencil[0].pos
            - stencil[0].pos*stencil[0].pos*stencil[1].pos
            - stencil[0].pos*stencil[0].pos*stencil[2].pos
            - stencil[0].pos*stencil[0].pos*stencil[3].pos
            + stencil[0].pos*stencil[1].pos*stencil[2].pos
            + stencil[0].pos*stencil[1].pos*stencil[3].pos
            + stencil[0].pos*stencil[2].pos*stencil[3].pos
            - stencil[1].pos*stencil[2].pos*stencil[3].pos);

  real c1 = -(stencil[0].pos*stencil[2].pos
            + stencil[0].pos*stencil[3].pos
            + stencil[2].pos*stencil[3].pos)
            /(stencil[0].pos*stencil[1].pos*stencil[1].pos
            - stencil[0].pos*stencil[1].pos*stencil[2].pos
            - stencil[0].pos*stencil[1].pos*stencil[3].pos
            + stencil[0].pos*stencil[2].pos*stencil[3].pos
            - stencil[1].pos*stencil[1].pos*stencil[1].pos
            + stencil[1].pos*stencil[1].pos*stencil[2].pos
            + stencil[1].pos*stencil[1].pos*stencil[3].pos
            - stencil[1].pos*stencil[2].pos*stencil[3].pos);

  real c2 =  (stencil[0].pos*stencil[1].pos
            + stencil[0].pos*stencil[3].pos
            + stencil[1].pos*stencil[3].pos)
            /(stencil[0].pos*stencil[1].pos*stencil[2].pos
            - stencil[0].pos*stencil[1].pos*stencil[3].pos
            - stencil[0].pos*stencil[2].pos*stencil[2].pos
            + stencil[0].pos*stencil[2].pos*stencil[3].pos
            - stencil[1].pos*stencil[2].pos*stencil[2].pos
            + stencil[1].pos*stencil[2].pos*stencil[3].pos
            + stencil[2].pos*stencil[2].pos*stencil[2].pos
            - stencil[2].pos*stencil[2].pos*stencil[3].pos);

  real c3 = -(stencil[0].pos*stencil[1].pos
            + stencil[0].pos*stencil[2].pos
            + stencil[1].pos*stencil[2].pos)
            /(stencil[0].pos*stencil[1].pos*stencil[2].pos
            - stencil[0].pos*stencil[1].pos*stencil[3].pos
            - stencil[0].pos*stencil[2].pos*stencil[3].pos
            + stencil[0].pos*stencil[3].pos*stencil[3].pos
            - stencil[1].pos*stencil[2].pos*stencil[3].pos
            + stencil[1].pos*stencil[3].pos*stencil[3].pos
            + stencil[2].pos*stencil[3].pos*stencil[3].pos
            - stencil[3].pos*stencil[3].pos*stencil[3].pos);
  
  coefs = { c0, c1, c2, c3 };
  return;
}

/******************************************************************************/
void Topology::fourth_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a fourth-order difference.
*******************************************************************************/
  real c0 = -(stencil[1].pos*stencil[2].pos*stencil[3].pos
            + stencil[1].pos*stencil[2].pos*stencil[4].pos
            + stencil[1].pos*stencil[3].pos*stencil[4].pos
            + stencil[2].pos*stencil[3].pos*stencil[4].pos)
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

  real c1 =  (stencil[0].pos*stencil[2].pos*stencil[3].pos
            + stencil[0].pos*stencil[2].pos*stencil[4].pos
            + stencil[0].pos*stencil[3].pos*stencil[4].pos
            + stencil[2].pos*stencil[3].pos*stencil[4].pos)
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

  real c2 = -(stencil[0].pos*stencil[1].pos*stencil[3].pos
            + stencil[0].pos*stencil[1].pos*stencil[4].pos
            + stencil[0].pos*stencil[3].pos*stencil[4].pos
            + stencil[1].pos*stencil[3].pos*stencil[4].pos)
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

  real c3 =  (stencil[0].pos*stencil[1].pos*stencil[2].pos
            + stencil[0].pos*stencil[1].pos*stencil[4].pos
            + stencil[0].pos*stencil[2].pos*stencil[4].pos
            + stencil[1].pos*stencil[2].pos*stencil[4].pos)
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

  real c4 = -(stencil[0].pos*stencil[1].pos*stencil[2].pos
            + stencil[0].pos*stencil[1].pos*stencil[3].pos
            + stencil[0].pos*stencil[2].pos*stencil[3].pos
            + stencil[1].pos*stencil[2].pos*stencil[3].pos)
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
  
  coefs = { c0, c1, c2, c3, c4 };
  return;
}

#else /* the functions below assume that zeroth point has stencil[0].pos = 0 */ 
/******************************************************************************/
void Topology::second_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1].pos - 1./stencil[2].pos;
  real c1 = -stencil[2].pos/(stencil[1].pos*(stencil[1].pos-stencil[2].pos));
  real c2 =  stencil[1].pos/(stencil[2].pos*(stencil[1].pos-stencil[2].pos));

  coefs = { c0, c1, c2 };
  return;
}

/******************************************************************************/
void Topology::third_order_first_coefs(std::vector<real> & coefs,
                                 const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a third-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1].pos - 1./stencil[2].pos - 1./stencil[3].pos; 
  
  real c1 =  stencil[2].pos*stencil[3].pos
           /(stencil[1].pos*(stencil[1].pos*stencil[1].pos
                       - stencil[1].pos*stencil[2].pos
                       - stencil[1].pos*stencil[3].pos
                       + stencil[2].pos*stencil[3].pos));

  real c2 = -stencil[1].pos*stencil[3].pos
           /(stencil[2].pos*(stencil[1].pos*stencil[2].pos
                       - stencil[1].pos*stencil[3].pos
                       - stencil[2].pos*stencil[2].pos
                       + stencil[2].pos*stencil[3].pos));

  real c3 =  stencil[1].pos*stencil[2].pos
           /(stencil[3].pos*(stencil[1].pos*stencil[2].pos
                       - stencil[1].pos*stencil[3].pos
                       - stencil[2].pos*stencil[3].pos
                       + stencil[3].pos*stencil[3].pos));

  coefs = { c0, c1, c2, c3 };
  return;
}

/******************************************************************************/
void Topology::fourth_order_first_coefs(std::vector<real> & coefs,
                                  const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a fourth-order difference.
*******************************************************************************/
  
  real c0 = -1./stencil[1].pos - 1./stencil[2].pos
           - 1./stencil[3].pos - 1./stencil[4].pos; 

  real c1 = -stencil[2].pos*stencil[3].pos*stencil[4].pos
           /(stencil[1].pos*(stencil[1].pos*stencil[1].pos*stencil[1].pos
                       - stencil[1].pos*stencil[1].pos*stencil[2].pos
                       - stencil[1].pos*stencil[1].pos*stencil[3].pos
                       - stencil[1].pos*stencil[1].pos*stencil[4].pos
                       + stencil[1].pos*stencil[2].pos*stencil[3].pos
                       + stencil[1].pos*stencil[2].pos*stencil[4].pos
                       + stencil[1].pos*stencil[3].pos*stencil[4].pos
                       - stencil[2].pos*stencil[3].pos*stencil[4].pos));

  real c2 =  stencil[1].pos*stencil[3].pos*stencil[4].pos
           /(stencil[2].pos*(stencil[1].pos*stencil[2].pos*stencil[2].pos
                       - stencil[1].pos*stencil[2].pos*stencil[3].pos
                       - stencil[1].pos*stencil[2].pos*stencil[4].pos
                       + stencil[1].pos*stencil[3].pos*stencil[4].pos
                       - stencil[2].pos*stencil[2].pos*stencil[2].pos
                       + stencil[2].pos*stencil[2].pos*stencil[3].pos
                       + stencil[2].pos*stencil[2].pos*stencil[4].pos
                       - stencil[2].pos*stencil[3].pos*stencil[4].pos));

  real c3 = -stencil[1].pos*stencil[2].pos*stencil[4].pos
           /(stencil[3].pos*(stencil[1].pos*stencil[2].pos*stencil[3].pos
                       - stencil[1].pos*stencil[2].pos*stencil[4].pos
                       - stencil[1].pos*stencil[3].pos*stencil[3].pos
                       + stencil[1].pos*stencil[3].pos*stencil[4].pos
                       - stencil[2].pos*stencil[3].pos*stencil[3].pos
                       + stencil[2].pos*stencil[3].pos*stencil[4].pos
                       + stencil[3].pos*stencil[3].pos*stencil[3].pos
                       - stencil[3].pos*stencil[3].pos*stencil[4].pos));

  real c4 =  stencil[1].pos*stencil[2].pos*stencil[3].pos
           /(stencil[4].pos*(stencil[1].pos*stencil[2].pos*stencil[3].pos
                       - stencil[1].pos*stencil[2].pos*stencil[4].pos
                       - stencil[1].pos*stencil[3].pos*stencil[4].pos
                       + stencil[1].pos*stencil[4].pos*stencil[4].pos
                       - stencil[2].pos*stencil[3].pos*stencil[4].pos
                       + stencil[2].pos*stencil[4].pos*stencil[4].pos
                       + stencil[3].pos*stencil[4].pos*stencil[4].pos
                       - stencil[4].pos*stencil[4].pos*stencil[4].pos));
  
  coefs = { c0, c1, c2, c3, c4 };
  return;
}
#endif
