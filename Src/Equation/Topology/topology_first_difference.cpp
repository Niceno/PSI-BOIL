#include "topology.h"

/* 
   presented expressions were calculated using symbolic python package sympy,
   based on theory from 
     https://en.wikipedia.org/wiki/Finite_difference_coefficient#cite_note-4
     http://web.media.mit.edu/~crtaylor/calculator.html
*/

/******************************************************************************/
real Topology::nth_order_first(const std::vector<real> & stencil,
                               const std::vector<real> & values,
                               const AccuracyOrder & order) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a desired-order difference.
*******************************************************************************/
  switch(order.eval()) {
    case 0 :
      return zeroth_order_first(stencil,values);
    case 1 :
      return first_order_first(stencil,values);
    case 2 :
      return second_order_first(stencil,values);
    case 3 :
      return third_order_first(stencil,values);
    case 4 :
      return fourth_order_first(stencil,values);
    default :
      boil::aout<<"Topology: unrecognised first difference requested. Exiting."
                <<boil::endl;
      exit(0);
  }

  return 0.0;
}

/******************************************************************************/
real Topology::zeroth_order_first(const std::vector<real> & stencil,
                                  const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a zeroth-order difference.
*******************************************************************************/
  return 0.0;
}

/******************************************************************************/
real Topology::first_order_first(const std::vector<real> & stencil,
                                 const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a first-order difference.
*******************************************************************************/
  real c0 =  1./(stencil[0]-stencil[1]);

  real c1 = -1./(stencil[0]-stencil[1]);

  return c0*values[0]+c1*values[1];
}

#if 1
/******************************************************************************/
real Topology::second_order_first(const std::vector<real> & stencil,
                                  const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  real c0 = -(stencil[1] + stencil[2])
            /(stencil[0]*stencil[0]
            - stencil[0]*stencil[1]
            - stencil[0]*stencil[2]
            + stencil[1]*stencil[2]);

  real c1 =  (stencil[0] + stencil[2])
            /(stencil[0]*stencil[1]
            - stencil[0]*stencil[2]
            - stencil[1]*stencil[1]
            + stencil[1]*stencil[2]);

  real c2 = -(stencil[0] + stencil[1])
            /(stencil[0]*stencil[1]
            - stencil[0]*stencil[2]
            - stencil[1]*stencil[2]
            + stencil[2]*stencil[2]);

  return c0*values[0]+c1*values[1]+c2*values[2];
}

/******************************************************************************/
real Topology::third_order_first(const std::vector<real> & stencil,
                                 const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a third-order difference.
*******************************************************************************/
  real c0 =  (stencil[1]*stencil[2] 
                       + stencil[1]*stencil[3] + stencil[2]*stencil[3])
            /(stencil[0]*stencil[0]*stencil[0]
            - stencil[0]*stencil[0]*stencil[1]
            - stencil[0]*stencil[0]*stencil[2]
            - stencil[0]*stencil[0]*stencil[3]
            + stencil[0]*stencil[1]*stencil[2]
            + stencil[0]*stencil[1]*stencil[3]
            + stencil[0]*stencil[2]*stencil[3]
            - stencil[1]*stencil[2]*stencil[3]);

  real c1 = -(stencil[0]*stencil[2]
                       + stencil[0]*stencil[3] + stencil[2]*stencil[3])
            /(stencil[0]*stencil[1]*stencil[1]
            - stencil[0]*stencil[1]*stencil[2]
            - stencil[0]*stencil[1]*stencil[3]
            + stencil[0]*stencil[2]*stencil[3]
            - stencil[1]*stencil[1]*stencil[1]
            + stencil[1]*stencil[1]*stencil[2]
            + stencil[1]*stencil[1]*stencil[3]
            - stencil[1]*stencil[2]*stencil[3]);

  real c2 =  (stencil[0]*stencil[1]
                       + stencil[0]*stencil[3] + stencil[1]*stencil[3])
            /(stencil[0]*stencil[1]*stencil[2]
            - stencil[0]*stencil[1]*stencil[3]
            - stencil[0]*stencil[2]*stencil[2]
            + stencil[0]*stencil[2]*stencil[3]
            - stencil[1]*stencil[2]*stencil[2]
            + stencil[1]*stencil[2]*stencil[3]
            + stencil[2]*stencil[2]*stencil[2]
            - stencil[2]*stencil[2]*stencil[3]);

  real c3 = -(stencil[0]*stencil[1]
                       + stencil[0]*stencil[2] + stencil[1]*stencil[2])
            /(stencil[0]*stencil[1]*stencil[2]
            - stencil[0]*stencil[1]*stencil[3]
            - stencil[0]*stencil[2]*stencil[3]
            + stencil[0]*stencil[3]*stencil[3]
            - stencil[1]*stencil[2]*stencil[3]
            + stencil[1]*stencil[3]*stencil[3]
            + stencil[2]*stencil[3]*stencil[3]
            - stencil[3]*stencil[3]*stencil[3]);
  
  return c0*values[0]+c1*values[1]+c2*values[2]+c3*values[3];
}

/******************************************************************************/
real Topology::fourth_order_first(const std::vector<real> & stencil,
                                  const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a fourth-order difference.
*******************************************************************************/
  real c0 = -(stencil[1]*stencil[2]*stencil[3] + stencil[1]*stencil[2]*stencil[4]
            + stencil[1]*stencil[3]*stencil[4] + stencil[2]*stencil[3]*stencil[4])
            /(stencil[0]*stencil[0]*stencil[0]*stencil[0]
            - stencil[0]*stencil[0]*stencil[0]*stencil[1]
            - stencil[0]*stencil[0]*stencil[0]*stencil[2]
            - stencil[0]*stencil[0]*stencil[0]*stencil[3]
            - stencil[0]*stencil[0]*stencil[0]*stencil[4]
            + stencil[0]*stencil[0]*stencil[1]*stencil[2]
            + stencil[0]*stencil[0]*stencil[1]*stencil[3]
            + stencil[0]*stencil[0]*stencil[1]*stencil[4]
            + stencil[0]*stencil[0]*stencil[2]*stencil[3]
            + stencil[0]*stencil[0]*stencil[2]*stencil[4]
            + stencil[0]*stencil[0]*stencil[3]*stencil[4]
            - stencil[0]*stencil[1]*stencil[2]*stencil[3]
            - stencil[0]*stencil[1]*stencil[2]*stencil[4]
            - stencil[0]*stencil[1]*stencil[3]*stencil[4] 
            - stencil[0]*stencil[2]*stencil[3]*stencil[4]
            + stencil[1]*stencil[2]*stencil[3]*stencil[4]);

  real c1 =  (stencil[0]*stencil[2]*stencil[3] + stencil[0]*stencil[2]*stencil[4]
            + stencil[0]*stencil[3]*stencil[4] + stencil[2]*stencil[3]*stencil[4])
            /(stencil[0]*stencil[1]*stencil[1]*stencil[1]
            - stencil[0]*stencil[1]*stencil[1]*stencil[2]
            - stencil[0]*stencil[1]*stencil[1]*stencil[3]
            - stencil[0]*stencil[1]*stencil[1]*stencil[4]
            + stencil[0]*stencil[1]*stencil[2]*stencil[3]
            + stencil[0]*stencil[1]*stencil[2]*stencil[4]
            + stencil[0]*stencil[1]*stencil[3]*stencil[4]
            - stencil[0]*stencil[2]*stencil[3]*stencil[4]
            - stencil[1]*stencil[1]*stencil[1]*stencil[1]
            + stencil[1]*stencil[1]*stencil[1]*stencil[2]
            + stencil[1]*stencil[1]*stencil[1]*stencil[3]
            + stencil[1]*stencil[1]*stencil[1]*stencil[4]
            - stencil[1]*stencil[1]*stencil[2]*stencil[3]
            - stencil[1]*stencil[1]*stencil[2]*stencil[4]
            - stencil[1]*stencil[1]*stencil[3]*stencil[4]
            + stencil[1]*stencil[2]*stencil[3]*stencil[4]);

  real c2 = -(stencil[0]*stencil[1]*stencil[3] + stencil[0]*stencil[1]*stencil[4]
            + stencil[0]*stencil[3]*stencil[4] + stencil[1]*stencil[3]*stencil[4])
            /(stencil[0]*stencil[1]*stencil[2]*stencil[2]
            - stencil[0]*stencil[1]*stencil[2]*stencil[3]
            - stencil[0]*stencil[1]*stencil[2]*stencil[4]
            + stencil[0]*stencil[1]*stencil[3]*stencil[4]
            - stencil[0]*stencil[2]*stencil[2]*stencil[2]
            + stencil[0]*stencil[2]*stencil[2]*stencil[3]
            + stencil[0]*stencil[2]*stencil[2]*stencil[4]
            - stencil[0]*stencil[2]*stencil[3]*stencil[4]
            - stencil[1]*stencil[2]*stencil[2]*stencil[2]
            + stencil[1]*stencil[2]*stencil[2]*stencil[3]
            + stencil[1]*stencil[2]*stencil[2]*stencil[4]
            - stencil[1]*stencil[2]*stencil[3]*stencil[4]
            + stencil[2]*stencil[2]*stencil[2]*stencil[2]
            - stencil[2]*stencil[2]*stencil[2]*stencil[3]
            - stencil[2]*stencil[2]*stencil[2]*stencil[4]
            + stencil[2]*stencil[2]*stencil[3]*stencil[4]);

  real c3 =  (stencil[0]*stencil[1]*stencil[2] + stencil[0]*stencil[1]*stencil[4]
            + stencil[0]*stencil[2]*stencil[4] + stencil[1]*stencil[2]*stencil[4])
            /(stencil[0]*stencil[1]*stencil[2]*stencil[3]
            - stencil[0]*stencil[1]*stencil[2]*stencil[4]
            - stencil[0]*stencil[1]*stencil[3]*stencil[3]
            + stencil[0]*stencil[1]*stencil[3]*stencil[4]
            - stencil[0]*stencil[2]*stencil[3]*stencil[3]
            + stencil[0]*stencil[2]*stencil[3]*stencil[4]
            + stencil[0]*stencil[3]*stencil[3]*stencil[3]
            - stencil[0]*stencil[3]*stencil[3]*stencil[4]
            - stencil[1]*stencil[2]*stencil[3]*stencil[3]
            + stencil[1]*stencil[2]*stencil[3]*stencil[4]
            + stencil[1]*stencil[3]*stencil[3]*stencil[3]
            - stencil[1]*stencil[3]*stencil[3]*stencil[4]
            + stencil[2]*stencil[3]*stencil[3]*stencil[3]
            - stencil[2]*stencil[3]*stencil[3]*stencil[4]
            - stencil[3]*stencil[3]*stencil[3]*stencil[3]
            + stencil[3]*stencil[3]*stencil[3]*stencil[4]);

  real c4 = -(stencil[0]*stencil[1]*stencil[2] + stencil[0]*stencil[1]*stencil[3]
            + stencil[0]*stencil[2]*stencil[3] + stencil[1]*stencil[2]*stencil[3])
            /(stencil[0]*stencil[1]*stencil[2]*stencil[3]
            - stencil[0]*stencil[1]*stencil[2]*stencil[4]
            - stencil[0]*stencil[1]*stencil[3]*stencil[4]
            + stencil[0]*stencil[1]*stencil[4]*stencil[4]
            - stencil[0]*stencil[2]*stencil[3]*stencil[4]
            + stencil[0]*stencil[2]*stencil[4]*stencil[4]
            + stencil[0]*stencil[3]*stencil[4]*stencil[4]
            - stencil[0]*stencil[4]*stencil[4]*stencil[4]
            - stencil[1]*stencil[2]*stencil[3]*stencil[4]
            + stencil[1]*stencil[2]*stencil[4]*stencil[4]
            + stencil[1]*stencil[3]*stencil[4]*stencil[4]
            - stencil[1]*stencil[4]*stencil[4]*stencil[4]
            + stencil[2]*stencil[3]*stencil[4]*stencil[4]
            - stencil[2]*stencil[4]*stencil[4]*stencil[4]
            - stencil[3]*stencil[4]*stencil[4]*stencil[4]
            + stencil[4]*stencil[4]*stencil[4]*stencil[4]);
  
  return c0*values[0]+c1*values[1]+c2*values[2]+c3*values[3]+c4*values[4];
}

#else /* the functions below assume that zeroth point has stencil[0] = 0 */ 
/******************************************************************************/
real Topology::second_order_first(const std::vector<real> & stencil,
                                  const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a second-order difference.
*******************************************************************************/
  real c0 = -1./stencil[1] - 1./stencil[2];
  real c1 = -stencil[2]/(stencil[1]*(stencil[1]-stencil[2]));
  real c2 =  stencil[1]/(stencil[2]*(stencil[1]-stencil[2]));

  return c0*values[0]+c1*values[1]+c2*values[2];
}

/******************************************************************************/
real Topology::third_order_first(const std::vector<real> & stencil,
                                 const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a third-order difference.
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
real Topology::fourth_order_first(const std::vector<real> & stencil,
                                  const std::vector<real> & values) const {
/***************************************************************************//**
*  \brief Approximate first derivative using a fourth-order difference.
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
#endif
