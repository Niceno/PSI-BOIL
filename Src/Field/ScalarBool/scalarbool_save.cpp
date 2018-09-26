#include "scalarbool.h"

/******************************************************************************/
void ScalarBool::save(const char * nm, const int it) {

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
void ScalarBool::save(std::ofstream & out) {

  int n_x = ni();
  int n_y = nj();
  int n_z = nk();
  
  out.write(reinterpret_cast<const char *> (&n_x), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_y), sizeof(int));
  out.write(reinterpret_cast<const char *> (&n_z), sizeof(int));

  out.write(reinterpret_cast<const char *> (val[0][0]), 
            ni()*nj()*nk()*sizeof(bool));
}

/*-----------------------------------------------------------------------------+
 '$Id: scalarbool_save.cpp,v 1.1 2014/02/04 08:16:57 sato Exp $'/
+-----------------------------------------------------------------------------*/
