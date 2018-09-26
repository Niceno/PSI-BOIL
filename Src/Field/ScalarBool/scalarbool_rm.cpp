#include "scalarbool.h"

/******************************************************************************/
void ScalarBool::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarbool_rm.cpp,v 1.1 2014/02/04 08:16:57 sato Exp $'/
+-----------------------------------------------------------------------------*/
