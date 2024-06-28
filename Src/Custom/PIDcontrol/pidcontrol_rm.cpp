#include "pidcontrol.h"

/***************************************************************************//**
*  save PIDcontrol data
*******************************************************************************/
void PIDcontrol::rm(const char * nm, const int it){
  if (boil::cart.iam()==0) {
    /* file name */
    std::string name = name_file(nm, ".bck", it);

    /* remove a file */
    remove(name.c_str());
  }
}

