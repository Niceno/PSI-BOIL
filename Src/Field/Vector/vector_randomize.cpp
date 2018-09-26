#include "vector.h"                

/******************************************************************************/
void Vector::randomize(RandomFlow & rf, const real magn, const real time) {

  boil::timer.start("randomization");

  for_m(m)
    for_amijk(m,i,j,k)
      vec[m][i][j][k]=0.0;

  /*-------------------------------+
  |  create random velocity field  |
  +-------------------------------*/
  real max=0.0;
  real x[3];
  for_m(m) {

/* theoretically, it should be like this */
//  for(int i=si(m); i<=ei(m); i++)
//  for(int j=sj(m); j<=ej(m); j++)
//  for(int k=sk(m); k<=ek(m); k++) {

    int sii = si(m);
    int eii = ei(m);
    int sjj = sj(m);
    int ejj = ej(m);
    int skk = sk(m);
    int ekk = ek(m);
    if(vec[m].bc().type_here(Dir::imin(),BndType::periodic())) sii=si(m)+3;
    if(vec[m].bc().type_here(Dir::imax(),BndType::periodic())) eii=ei(m)-3;
    if(vec[m].bc().type_here(Dir::jmin(),BndType::periodic())) sjj=sj(m)+3;
    if(vec[m].bc().type_here(Dir::jmax(),BndType::periodic())) ejj=ej(m)-3;
    if(vec[m].bc().type_here(Dir::kmin(),BndType::periodic())) skk=sk(m)+3;
    if(vec[m].bc().type_here(Dir::kmax(),BndType::periodic())) ekk=ek(m)-3;

    /* crude code */
    for(int i=sii; i<=eii; i++)
    for(int j=sjj; j<=ejj; j++)
    for(int k=skk; k<=ekk; k++) {
      x[0] = xc(m,i); 
      x[1] = yc(m,j); 
      x[2] = zc(m,k);
      rf.get_vector( time, x, &vec[m][i][j][k], m );
      max = boil::maxr(fabs(vec[m][i][j][k]), max);
    }
  }

  /*-------------------------+ 
  |  scale the fluctuations  |
  +-------------------------*/
  real factor = magn/max;
  for_m(m) {
    for_mijk(m,i,j,k) {
      vec[m][i][j][k] *= factor;
    }
  }

  /*-----------------------------------------------------+
  |  distribute velocities from first to all processors  |
  +-----------------------------------------------------*/
#if 0
  vec[Comp::u()].one2all(0);  
  vec[Comp::v()].one2all(0);  
  vec[Comp::w()].one2all(0);  
#endif

  boil::timer.stop("randomization");
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_randomize.cpp,v 1.15 2015/09/23 09:24:40 sato Exp $'/
+-----------------------------------------------------------------------------*/
