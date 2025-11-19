#include "vector.h"

/******************************************************************************/
void Vector::save(const char * nm, const int it) {

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
void Vector::save(std::ofstream & out) {

  for_m(m)
   {int n_x = ni(m);
    int n_y = nj(m);
    int n_z = nk(m);
    out.write(reinterpret_cast<const char *> (&n_x), sizeof(int));
    out.write(reinterpret_cast<const char *> (&n_y), sizeof(int));
    out.write(reinterpret_cast<const char *> (&n_z), sizeof(int));}

  for_m(m)
    out.write(reinterpret_cast<const char *> (vec[m].val[0][0]), 
              ni(m)*nj(m)*nk(m) * sizeof(real));
}
