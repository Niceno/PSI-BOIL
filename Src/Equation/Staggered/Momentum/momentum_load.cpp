#include "momentum.h"

/******************************************************************************/
void Momentum::load(const char * nm, const int it) {

  std::string name;

  /* file name */
  if(it == -1) 
    name = name_file(nm, ".bck", time->first_step(), boil::cart.iam());
  else
    name = name_file(nm, ".bck", it,                 boil::cart.iam());

  /* open a file */
  std::ifstream in(name.c_str(), std::ios::binary);
  
  /* stop if file is not present */
  if( in.rdstate() != 0 ) {
    std::cout << "failed to open " << name << std::endl;
    std::cout << "exiting!" << std::endl;
    exit(0);
  }

  /* load the necessary data */
  u.   load(in);
  cold.load(in);

  /* close a file */
  in.close();
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_load.cpp,v 1.5 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
