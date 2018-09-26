#include "vectorbool.h"

/******************************************************************************/
void VectorBool::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* remove file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: vectorbool_rm.cpp,v 1.1 2014/02/04 08:20:58 sato Exp $'/
+-----------------------------------------------------------------------------*/
