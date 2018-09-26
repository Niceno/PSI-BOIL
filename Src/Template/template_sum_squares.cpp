#include "template.h" /* class header must be included */

/***************************************************************************//**
*  This function mimics a definition of a non-trivial member function.
*******************************************************************************/
real Template :: sum_squares() const {

  /*-----------------------------------+
  |  compute the result and return it  |
  +-----------------------------------*/
  return at*at + bt*bt;
}

/*-----------------------------------------------------------------------------+
 '$Id: template_sum_squares.cpp,v 1.1 2009/11/12 15:04:28 niceno Exp $'/
+-----------------------------------------------------------------------------*/
