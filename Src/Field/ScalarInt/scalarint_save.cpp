#include "scalarint.h"

/******************************************************************************/
void ScalarInt::save(const char * nm, const int it) {

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
void ScalarInt::save(std::ofstream & out) {

  int n_x = ni();
  int n_y = nj();
  int n_z = nk();
  
  out.write(reinterpret_cast<const char *> (&n_x), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_y), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_z), sizeof(int));

  out.write(reinterpret_cast<const char *> (val[0][0]), 
            ni()*nj()*nk()*sizeof(int));
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarint_save.cpp,v 1.1 2015/05/05 14:36:01 sato Exp $'/
+-----------------------------------------------------------------------------*/
