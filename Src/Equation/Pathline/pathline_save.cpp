#include "pathline.h"

/***************************************************************************//**
*  save pathline data
*******************************************************************************/
void Pathline::save(const char * nm, const int it) {

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
void Pathline::save(std::ofstream & out) {

  int nparticles = np();
  out.write(reinterpret_cast<const char *> (&nparticles), sizeof(int));

  for (int ip = 0; ip < np(); ip++){
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    int  idd  = particles[ip].id();
    out.write(reinterpret_cast<const char *> (&xold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&yold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&zold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&idd ), sizeof(int));
  }
  return;
}


