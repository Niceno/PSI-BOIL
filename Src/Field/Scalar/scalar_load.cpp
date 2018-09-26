#include "scalar.h"

/******************************************************************************/
void Scalar::load(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  std::ifstream in(name.c_str(), std::ios::binary);
  
  /* stop if file is not present */
  if( in.rdstate() != 0 ) {
    std::cout << "failed to open " << name << std::endl;
    std::cout << "exiting!" << std::endl;
    exit(0);
  }

  /* load the necessary data */
  load(in);
  
  /* close a file */
  in.close();
}

/******************************************************************************/
void Scalar::load(std::ifstream & in) {

  int n_x_saved, n_y_saved, n_z_saved;
  
  in.read(reinterpret_cast<char *> (&n_x_saved), sizeof(int));
  in.read(reinterpret_cast<char *> (&n_y_saved), sizeof(int));
  in.read(reinterpret_cast<char *> (&n_z_saved), sizeof(int));

  in.read(reinterpret_cast<char *> (val[0][0]), 
          ni()*nj()*nk()*sizeof(real));
}
/*-----------------------------------------------------------------------------+
 '$Id: scalar_load.cpp,v 1.9 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
