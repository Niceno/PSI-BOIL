#include "solidparticle.h"

/***************************************************************************//**
*  add solid particle
*******************************************************************************/
void SolidParticle::add_global(const real x, const real y, const real z,
                               const real dia, const real den) {

#if 0
  // check input
  real xx = x;
  real yy = y;
  real zz = z;
  boil::cart.sum_real(&xx);
  boil::cart.sum_real(&yy);
  boil::cart.sum_real(&zz);
  int nproc = boil::cart.nproc();
  if (!approx(x,xx/real(nproc),boil::pico) ||
      !approx(y,yy/real(nproc),boil::pico) ||
      !approx(z,zz/real(nproc),boil::pico)) {
    std::cout<<"pathline_add:ERROR!!!\n";
    std::cout<<"You need to use pathline_add_local and pathline.exchange\n";
    exit(0);
  }
#endif

  // define id
  id_serial++;

  // create particle
  //Particle p(x, y, z, nval());
  //Particle p(x, y, z, 1);
  Particle p(x, y, z, nval(), dia, den);
  //Particle p(x, y, z, dia, den, 0);

  // give id to particle
  p.id(id_serial);

  // push_back
  particles.push_back(p);

  // get np
  np(particles.size());

  if (np()%10==1)
  boil::oout<<"Pathline:add:np= "<<npa<<" x= "<<particles[npa-1].x()<<" y= "
            <<particles[npa-1].y()<<" z= "<<particles[npa-1].z()
            <<" dia= "<<particles[npa-1].diameter()
            <<" den= "<<particles[npa-1].density()<<"\n";
}

