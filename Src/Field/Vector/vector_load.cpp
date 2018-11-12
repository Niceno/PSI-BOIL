#include "vector.h"

/******************************************************************************/
void Vector::load(const char * nm, const int it) {

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
void Vector::load(std::ifstream & in) {

  int n_x_saved[3], n_y_saved[3], n_z_saved[3];
  
  for_m(m)
   {in.read(reinterpret_cast<char *> (&n_x_saved[~m]), sizeof(int));
    in.read(reinterpret_cast<char *> (&n_y_saved[~m]), sizeof(int));
    in.read(reinterpret_cast<char *> (&n_z_saved[~m]), sizeof(int));}

  for_m(m)
    in.read(reinterpret_cast<char *> (vec[m].val[0][0]), 
            ni(m)*nj(m)*nk(m) * sizeof(real));
}
