#include "template.h" /* class header must be included */

/***************************************************************************//**
*  This is a non-trivial constructor. Note that this comment is written in
*  Doxygen style and therefore in a grammatically correct English language.            
*******************************************************************************/
Template :: Template(const real & x, const real & y) {

  /*----------------------------------------------+
  |  copy values of paremeters into data members  |
  +----------------------------------------------*/
  at = x;
  bt = y;
}

/******************************************************************************/
Template :: ~Template() {
/*--------------------------------------------------------------------------+
|  destructor is also placed here.                                          |
|  this is a psi-boil-style comment, written with lower-case letters only.  |
+--------------------------------------------------------------------------*/

}
