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
