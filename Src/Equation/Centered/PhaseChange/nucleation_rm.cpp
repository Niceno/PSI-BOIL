#include "nucleation.h"

/******************************************************************************/
void Nucleation::rm(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  remove(name.c_str());
  
}

/*-----------------------------------------------------------------------------+
 '$Id: nucleation_rm.cpp,v 1.1 2014/08/06 08:19:36 sato Exp $'/
+-----------------------------------------------------------------------------*/
