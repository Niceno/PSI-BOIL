#include "solidparticle.h"

/***************************************************************************//**
*  save solid particle data
*******************************************************************************/
void SolidParticle::load(const char * nm, const int it) {

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
void SolidParticle::load(std::ifstream & in) {

  int nparticles;
  real x_saved, y_saved, z_saved, u_saved, v_saved, w_saved, dia_saved, den_saved;
  int  i_saved;
  in.read(reinterpret_cast<char *> (&nparticles), sizeof(int));
  boil::oout<<"solidparticle_load: total number of particles = "<<nparticles<<"\n";

  for (int ip = 0; ip < nparticles; ip++){
    in.read(reinterpret_cast<char *> (&x_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&y_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&z_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&u_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&v_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&w_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&dia_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&den_saved), sizeof(real));
    in.read(reinterpret_cast<char *> (&i_saved), sizeof(int));
    add_global(x_saved, y_saved, z_saved);
    particles[ip].u(u_saved);
    particles[ip].v(v_saved);
    particles[ip].w(w_saved);
    particles[ip].diameter(dia_saved);
    particles[ip].density (den_saved);
    particles[ip].id(i_saved);
  }
  return;
}


