#include "scalar.h"

/******************************************************************************/
void Scalar::save(const char * nm, const int it) {

  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());

  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);
  
  /* save the necessary data */
  save(out);
  
  /* close a file */
  out.close();
}

/******************************************************************************/
void Scalar::save(std::ofstream & out) {

  int n_x = ni();
  int n_y = nj();
  int n_z = nk();
  
  out.write(reinterpret_cast<const char *> (&n_x), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_y), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_z), sizeof(int));

  out.write(reinterpret_cast<const char *> (val[0][0]), 
            ni()*nj()*nk()*sizeof(real));
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_save.cpp,v 1.8 2008/11/17 19:23:23 niceno Exp $'/
+-----------------------------------------------------------------------------*/
