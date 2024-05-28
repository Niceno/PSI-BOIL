#include "solidparticle.h"

/***************************************************************************//**
*  cynchronize particle_local
*******************************************************************************/
void SolidParticle::exchange() {

  // loop for processes
  for (int iproc=0; iproc<boil::cart.nproc(); iproc++) {
    // number of particles defined in decomposed domain
    int np_local=0; 
    if(boil::cart.iam()==iproc) {
      np_local = particle_local.size();
    }
    // number of local particles
    boil::cart.sum_int(&np_local);
    if (np_local !=0 ) {
      boil::oout<<"exchange:np_local= "<<iproc<<" "<<boil::cart.iam()<<" "
      <<np_local<<"\n";
    }

    // loop for local particles in iproc
    for (int ip=0; ip<np_local; ip++) {
      real xx=0.0, yy=0.0, zz=0.0, dia=0.0, den=0.0;
      if(boil::cart.iam()==iproc) {
        xx = particle_local[ip].x;
        yy = particle_local[ip].y;
        zz = particle_local[ip].z;
        dia = particle_local[ip].diameter;
        den = particle_local[ip].density;
      }
      boil::cart.sum_real(&xx);
      boil::cart.sum_real(&yy);
      boil::cart.sum_real(&zz);
      add_global(xx, yy, zz, dia, den);
      //boil::oout<<"add(xx, yy, zz): "<<xx<<" "<<yy<<" "<<zz<<"\n";
    }
  }
}

