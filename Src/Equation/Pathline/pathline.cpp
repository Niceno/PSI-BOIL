#include "pathline.h"

/***************************************************************************//**
*  Constructor for pathline
*******************************************************************************/
Pathline::Pathline ( const Vector & v, const Times * t,
                     const Scalar * sca1,
		     const Scalar * sca2,
		     const Scalar * sca3) {
  uvw = & v;
  time = t;
  s1 = sca1;
  s2 = sca2;
  s3 = sca3;
  npa = 0;
  nv = 0;
  id_serial = -1;
  if (s1 !=NULL) { nv = 1;}
  if (s2 !=NULL) { nv = 2;}
  if (s3 !=NULL) { nv = 3;}
  //std::cout<<"pathline:time= "<<time->current_time()<<"\n";
  boil::oout<<"Pathline:nval()= "<<nval()<<"\n";
}

/******************************************************************************/
Pathline::~Pathline() {
}
