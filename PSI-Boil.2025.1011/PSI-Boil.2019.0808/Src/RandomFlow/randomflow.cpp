#include "randomflow.h"

/******************************************************************************/
RandomFlow::RandomFlow(int Ne, real turb_time, real turb_length ) {
/*-------------------------------------------+
|  Generate spectral expansion coefficients  |
+-------------------------------------------*/
  ne=Ne;
  if (ne<=0) return;
  real fe=sqrt(2./(real)ne);

  alloc1d(&Omega,ne);
  alloc2d(&K ,ne,DIM);
  alloc2d(&U1,ne,DIM);
  alloc2d(&U2,ne,DIM);

  OMS(skipping seed to get the same sequence every time);
  boil::random_seed(1);

  gaussn(K,0.5,ne,DIM);

  for (int ie=0; ie<ne; ie++) {

    int i,j=ie*DIM;
    real a,V1[DIM],V2[DIM]; /* random vectors xi and zeta (Eq.16) */
    Omega[ie] = gauss()/turb_time; 
    gaussn(V1, fe, DIM);
    gaussn(V2, fe, DIM);

    /* U1 = V1 x K */
    /* U2 = V2 x K */
    vecp(U1[ie],V1,K[ie]); /* Eq.16\cite{LiAhetalJAS94} */
    vecp(U2[ie],V2,K[ie]);
    for(i=0; i<DIM; i++) 
      K[ie][i]/=turb_length;
  }
}

/******************************************************************************/
void RandomFlow::get_vector(const real t, const real *x,  real *v) {

  real a,c,s;

  for (int i=0; i<DIM; i++) 
    v[i]=0.0;

  if (ne<=0)  return;

  for (int ie=0; ie<ne; ie++) { 
    a=sclp(K[ie],x)+Omega[ie]*t;
    c=cos(a); s=sin(a);
    for (int i=0; i<DIM; i++) 
      v[i]+=U1[ie][i]*c+U2[ie][i]*s;
  }
}

/******************************************************************************/
void RandomFlow::get_vector(const real t, const real *x,  real *v, 
                            const Comp & i) { // component

  real a,c,s;

  *v=0.0;

  if (ne<=0)  return;

  for (int ie=0; ie<ne; ie++) { 
    a=sclp(K[ie],x)+Omega[ie]*t;
    c=cos(a); s=sin(a);
    *v+=U1[ie][~i]*c+U2[ie][~i]*s;
  }
}
