#include "centered.h"

/******************************************************************************/
void Centered::save(const char * nm, const int it) {

  std::string name;

  /* file name */
  if(it == -1) // "it" is not specified, use the time step
    name = name_file(nm, ".bck", time->current_step(), boil::cart.iam());
  else
    name = name_file(nm, ".bck", it,                   boil::cart.iam());

  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);
  
  /* save necessary variables */
  phi. save(out);
  cold.save(out);

  /* close a file */
  out.close();
}
