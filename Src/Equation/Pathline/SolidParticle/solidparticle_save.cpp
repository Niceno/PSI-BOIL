#include "solidparticle.h"

/***************************************************************************//**
*  save solid particle data
*******************************************************************************/
void SolidParticle::save(const char * nm, const int it) {

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
void SolidParticle::save(std::ofstream & out) {

  int nparticles = np();
  out.write(reinterpret_cast<const char *> (&nparticles), sizeof(int));

  for (int ip = 0; ip < np(); ip++){
    real xold = particles[ip].x();
    real yold = particles[ip].y();
    real zold = particles[ip].z();
    real uold = particles[ip].u();
    real vold = particles[ip].v();
    real wold = particles[ip].w();
    real dia  = particles[ip].diameter();
    real den  = particles[ip].density();
    int  idd  = particles[ip].id();
    out.write(reinterpret_cast<const char *> (&xold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&yold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&zold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&uold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&vold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&wold), sizeof(real));
    out.write(reinterpret_cast<const char *> (&dia), sizeof(real));
    out.write(reinterpret_cast<const char *> (&den), sizeof(real));
    out.write(reinterpret_cast<const char *> (&idd ), sizeof(int));
  }
  return;
}


