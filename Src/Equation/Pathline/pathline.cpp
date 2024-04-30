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

  boil::oout<<"###########################################################\n";
  boil::oout<<"#  Pathline: WARNING with respect to parallelization !!!  #\n";
  boil::oout<<"#  pathline_add(x, y, z) must be cynchronized between     #\n";
  boil::oout<<"#  processes. Next loop does not work                     #\n";
  boil::oout<<"#---------------------------------------------------------#\n";
  boil::oout<<"#    for_vijk(c,i,j,k) {                                  #\n";
  boil::oout<<"#        if(c[i][j][k]<0.5) {                             #\n";
  boil::oout<<"#          pathline.add_global(c.xc(i),c.yc(j),c.zc(k));  #\n";
  boil::oout<<"#        }                                                #\n";
  boil::oout<<"#    }                                                    #\n";
  boil::oout<<"###########################################################\n";
  boil::oout<<"#  Next is correct.                                       #\n";
  boil::oout<<"#---------------------------------------------------------#\n";
  boil::oout<<"#    for_vijk(c,i,j,k) {                                  #\n";
  boil::oout<<"#      if(c[i][j][k]<0.5) {                               #\n";
  boil::oout<<"#        pathline.add_local(c.xc(i),c.yc(j),c.zc(k));     #\n";
  boil::oout<<"#      }                                                  #\n";
  boil::oout<<"#    }                                                    #\n";
  boil::oout<<"#    pathline.exchange();                                 #\n";
  boil::oout<<"###########################################################\n";
  boil::oout<<"#  Next is also correct.                                  #\n";
  boil::oout<<"#---------------------------------------------------------#\n";
  boil::oout<<"#    pathline.add_global(0.001,0.001,0.001);              #\n";
  boil::oout<<"###########################################################\n";


}

/******************************************************************************/
Pathline::~Pathline() {
}
