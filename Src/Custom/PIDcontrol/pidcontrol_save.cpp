#include "pidcontrol.h"

/***************************************************************************//**
*  save PIDcontrol data
*******************************************************************************/
void PIDcontrol::save(const char * nm, const int it){
  if (boil::cart.iam()==0) {
    /* file name */
    //std::string name = name_file(nm, ".bck", it, boil::cart.iam());
    std::string name = name_file(nm, ".bck", it);

    /* open a file */
    std::ofstream out(name.c_str(), std::ios::binary);

    /* save the necessary data */
    save(out);

    /* close a file */
    out.close();
  }
}

/******************************************************************************/
void PIDcontrol::save(std::ofstream & out) {
      out.write(reinterpret_cast<const char *> (&last_time_), sizeof(real));
      out.write(reinterpret_cast<const char *> (&integral_), sizeof(real));
      out.write(reinterpret_cast<const char *> (&prev_error_), sizeof(real));
      return;
}

