#include "pathline.h"

/***************************************************************************//**
*  Constructor for pathline
*******************************************************************************/
Pathline::Pathline ( const Vector & v, const Times * t) {
  uvw = & v;
  time = t;
  npa = 0;
  //std::cout<<"pathline:time= "<<time->current_time()<<"\n";
}

/******************************************************************************/
Pathline::~Pathline() {
}
