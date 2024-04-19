#include "pathline.h"

/***************************************************************************//**
*  save pathline data
*******************************************************************************/
void Pathline::load(const char * nm, const int it) {

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
void Pathline::load(std::ifstream & in) {

  int nparticles;
  real x_saved, y_saved, z_saved;
  in.read(reinterpret_cast<char *> (&nparticles), sizeof(int));
  boil::oout<<"pathline_load: total number of particles = "<<nparticles<<"\n";

  for (int ip = 0; ip < nparticles; ip++){
    in.read(reinterpret_cast<char *> (&x_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&y_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&z_saved), sizeof(real));
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    add(x_saved, y_saved, z_saved);
  }
  return;
}


