#include "scalarint.h"

/******************************************************************************/
void ScalarInt::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarint_rm.cpp,v 1.1 2015/05/05 14:36:01 sato Exp $'/
+-----------------------------------------------------------------------------*/
