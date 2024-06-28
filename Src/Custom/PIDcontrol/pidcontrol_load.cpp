#include "pidcontrol.h"

/***************************************************************************//**
*  load PIDcontrol data
*******************************************************************************/
void PIDcontrol::load(const char * nm, const int it){

  // initialize
  last_time_  = 0.0;
  integral_   = 0.0;
  prev_error_ = 0.0;

  if (boil::cart.iam()==0) {
    /* file name */
    std::string name = name_file(nm, ".bck", it);

    /* open a file */
    std::ifstream in(name.c_str(), std::ios::binary);

    /* load the necessary data */
    load(in);

    /* close a file */
    in.close();
  }

  boil::cart.sum_real(&last_time_);
  boil::cart.sum_real(&integral_);
  boil::cart.sum_real(&prev_error_);

}

/******************************************************************************/
void PIDcontrol::load(std::ifstream & in) {
  in.read(reinterpret_cast<char *> (&last_time_), sizeof(real));
  in.read(reinterpret_cast<char *> (&integral_), sizeof(real));
  in.read(reinterpret_cast<char *> (&prev_error_), sizeof(real));
  return;
}

