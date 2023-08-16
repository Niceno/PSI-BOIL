#include "pathline.h"

/***************************************************************************//**
*  cynchronize particle_local
*******************************************************************************/
void Pathline::exchange() {
  //std::cout<<"exchange000:\n";

  // loop for processes
  for (int iproc=0; iproc<boil::cart.nproc(); iproc++) {
    // number of particles defined in decomposed domain
    int np_local=0; 
    if(boil::cart.iam()==iproc) {
      np_local = particle_local.size();
    }
    // number of local particles
    boil::cart.sum_int(&np_local);
    //if (np_local !=0 ) {
    //  boil::oout<<"exchange:np_local= "<<iproc<<" "<<boil::cart.iam()<<" "
    //  <<np_local<<"\n";
    //}

    // loop for local particles in iproc
    for (int ip=0; ip<np_local; ip++) {
      real xx=0.0, yy=0.0, zz=0.0, dia=0.0, den=0.0;
      if(boil::cart.iam()==iproc) {
        xx = particle_local[ip].x;
        yy = particle_local[ip].y;
        zz = particle_local[ip].z;
	//std::cout<<"xx "<<xx<<" "<<yy<<" "<<zz<<"\n";
	//std::cout<<"dia "<<particle_local[ip].diameter<<"\n";
	if (particle_local[ip].diameter!=NULL) {
          dia = *(particle_local[ip].diameter);
	}
	if (particle_local[ip].density!=NULL) {
          den = *(particle_local[ip].density);
	}
        //dia = *(particle_local[ip].diameter);
        //den = *(particle_local[ip].density);
	//std::cout<<"dia "<<dia<<" "<<den<<" "<<zz<<"\n";
      }
      boil::cart.sum_real(&xx);
      boil::cart.sum_real(&yy);
      boil::cart.sum_real(&zz);
      boil::cart.sum_real(&dia);
      boil::cart.sum_real(&den);
      //if(particle_local[0].diameter==NULL) {
      //std::cout<<"exchange110:\n";
      if(b_dia_den) {
        add_global(xx, yy, zz, &dia, &den);
	//std::cout<<"exchange:add_global(xx,yy,zz,dia,den)\n";
      } else {
        //std::cout<<"exchange120:\n";
        add_global(xx,yy,zz);
	//std::cout<<"exchange:add_global(xx,yy,zz)\n";
      }
      //boil::oout<<"add(xx, yy, zz): "<<xx<<" "<<yy<<" "<<zz<<"\n";
    }
  }
}

