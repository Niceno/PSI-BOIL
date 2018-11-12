#include "momentum.h"

/******************************************************************************/
void Momentum::save(const char * nm, const int it) {

  std::string name;

  /* file name */
  if(it == -1) 
    name = name_file(nm, ".bck", time->current_step(), boil::cart.iam());
  else
    name = name_file(nm, ".bck", it,                   boil::cart.iam());

  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);
  
  /* load the necessary data */
  u.   save(out);
  cold.save(out);

  /* close a file */
  out.close();
}
