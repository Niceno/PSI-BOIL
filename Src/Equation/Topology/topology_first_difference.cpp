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
  real c0 =  1./(stencil[0].pos-stencil[1].pos);

  real c1 = -1./(stencil[0].pos-stencil[1].pos);

  return c0*stencil[0].val+c1*stencil[1].val;
}

#if 1
/******************************************************************************/
real Topology::second_order_first(const std::vector<StencilPoint> & stencil)
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

  return c0*stencil[0].val+c1*stencil[1].val+c2*stencil[2].val;
}

/******************************************************************************/
real Topology::third_order_first(const std::vector<StencilPoint> & stencil)
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
  
  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val;
}

/******************************************************************************/
real Topology::fourth_order_first(const std::vector<StencilPoint> & stencil)
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
  
  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val+c4*stencil[4].val;
}

#else /* the functions below assume that zeroth point has stencil[0].pos = 0 */ 
/******************************************************************************/
real Topology::second_order_first(const std::vector<StencilPoint> & stencil)
                                                                         const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1].pos - 1./stencil[2].pos;
  real c1 = -stencil[2].pos/(stencil[1].pos*(stencil[1].pos-stencil[2].pos));
  real c2 =  stencil[1].pos/(stencil[2].pos*(stencil[1].pos-stencil[2].pos));

  return c0*stencil[0].val+c1*stencil[1].val+c2*stencil[2].val;
}

/******************************************************************************/
real Topology::third_order_first(const std::vector<StencilPoint> & stencil)
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

  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val;
}

/******************************************************************************/
real Topology::fourth_order_first(const std::vector<StencilPoint> & stencil)
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
  
  return c0*stencil[0].val+c1*stencil[1].val
        +c2*stencil[2].val+c3*stencil[3].val+c4*stencil[4].val;
}
#endif
