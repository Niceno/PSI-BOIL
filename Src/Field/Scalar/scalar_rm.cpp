#include "scalar.h"

/******************************************************************************/
void Scalar::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_rm.cpp,v 1.1 2011/09/22 14:56:05 sato Exp $'/
+-----------------------------------------------------------------------------*/
