// Compile: g++ hturb_main.cpp ../Global/global_malloc.cpp turb_common.cpp

#include "turb_common.h"

int ne = 1; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

real
  *Omega,    /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
  **U1,**U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
  **K;       /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */

/******************************************************************************/
void genspec(int *Ne, real *TURB_TIME, real *TURB_LENGTH ) {
/*-------------------------------------------+
|  Generate spectral expansion coefficients  |
+-------------------------------------------*/
  real fe, turb_time=*TURB_TIME;
  ne=*Ne;
  if (ne<=0) return;
  fe=sqrt(2./(real)ne);

  alloc1d(&Omega,ne);
  alloc2d(&K ,ne,DIM);
  alloc2d(&U1,ne,DIM);
  alloc2d(&U2,ne,DIM);

  seed(); /* initialize the random generator */

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
      K[ie][i]/=TURB_LENGTH[0];
  }
}

/******************************************************************************/
void genvec(const real t, const real *x,  real *v) {

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
 main() {
   real x[3],v[3];
   int nspec=1000;      // spectral sample size
   real turbtime=1.0;   // turbulence time-scale
   real turblength=1.0; // turbulence length-scale
 
   genspec(&nspec, &turbtime, &turblength);
 
   int  nt=100; // number of time steps
   real dt=0.5; // time-step
   int  nx= 32; // number space-points
   real dx=1.0; // their separation
   int  ny= 32; // number space-points
   real dy=1.0; // their separation
 
   real ** x2d;
   real ** y2d;
   real ** u2d;
   real ** v2d;
   real ** w2d;

   alloc2d(&x2d, nx, ny);
   alloc2d(&y2d, nx, ny);
   alloc2d(&u2d, nx, ny);
   alloc2d(&v2d, nx, ny);
   alloc2d(&w2d, nx, ny);

   for(int it=0; it<nt; it++) {

     /* file name extension */
     std::string        name("random");
     std::ostringstream numb;
     numb << ".gmv.";
     numb.fill('0');
     numb.width(4);
     numb << it;
     name += numb.str();

     /* open the file */
     FILE * fp = fopen(name.c_str(), "w");

     printf("Creating: %s\n", name.c_str());

     /* create and store random velocity fields */
     for(int ix=0; ix<nx; ix++) 
       for(int iy=0; iy<ny; iy++) {

         real t=dt*it;
         x[0]=dx*ix;
         x[1]=dy*iy;
         x[2]=0.0;
 
         genvec(t,x,v);

         x2d[ix][iy] = x[0];
         y2d[ix][iy] = x[1];

         u2d[ix][iy] = v[0];
         v2d[ix][iy] = v[1];
         w2d[ix][iy] = v[2];
       }

     fprintf(fp,"gmvinput ascii\n");
     fprintf(fp,"nodev %d\n", nx*ny);
     for(int ix=0; ix<nx; ix++) 
       for(int iy=0; iy<ny; iy++) 
         fprintf(fp, "%f %f %f\n", x2d[ix][iy], y2d[ix][iy], 0.0);

     fprintf(fp,"cells 0\n");
     fprintf(fp,"velocity 1\n");
     for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++) fprintf(fp, "%f\n", u2d[ix][iy]);
     for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++) fprintf(fp, "%f\n", v2d[ix][iy]);
     for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++) fprintf(fp, "%f\n", w2d[ix][iy]);

     fprintf(fp,"endgmv\n");

     fclose(fp);
   }
 }
