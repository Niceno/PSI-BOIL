// Compile: g++ nhturb_main.cpp ../Global/global_malloc.cpp turb_common.cpp jacobi.cpp

#include "turb_common.h"

void jacobi(real ** a, real * d, real ** v, const int n, int * nrot);

int ne = 1; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

real
  *Omega,    /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
  **U1,**U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
  **K;       /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */

/******************************************************************************/
void genspec(int *Ne) {
/*--------------------------------------------------------------------------+
|  Generate spectral expansion coefficients for non-homogeneous turbulence  |
+--------------------------------------------------------------------------*/
  ne=*Ne;
  if (ne<=0) return;

  alloc1d(&Omega,ne);
  alloc2d(&K ,ne,DIM);
  alloc2d(&U1,ne,DIM);
  alloc2d(&U2,ne,DIM);

//  seed(); /* initialize the random generator */

  gaussn(K,0.5,ne,DIM);

  for (int ie=0; ie<ne; ie++) { 

    int i,j=ie*DIM;
    real a,V1[DIM],V2[DIM]; /* random vectors xi and zeta (Eq.16) */
    Omega[ie] = gauss();  
    gaussn(V1, 1., DIM);
    gaussn(V2, 1., DIM);

    /* U1 = V1 x K */
    /* U2 = V2 x K */
    vecp(U1[ie],V1,K[ie]); /* Eq.16\cite{LiAhetalJAS94} */
    vecp(U2[ie],V2,K[ie]);
  }
}

/******************************************************************************/
void genvec(const real t, const real *x,  real tt, const real * uu_in, 
            real *v) {

  real a,c,s;
  real  k[DIM]; /* non-homogeneous k */
  real  ** uu, ** vv, * d;

  if (ne<=0)  return;

  /* allocate memory */
  alloc2d(&uu,DIM,DIM);
  alloc2d(&vv,DIM,DIM);
  alloc1d(&d, DIM);

  /* form uu tensor */
  uu[0][0] = uu_in[0]; // uu
  uu[1][1] = uu_in[1]; // vv
  uu[2][2] = uu_in[2]; // ww
  uu[0][1] = uu[1][0] = uu_in[3]; // uv
  uu[0][2] = uu[2][0] = uu_in[4]; // uw
  uu[1][2] = uu[2][1] = uu_in[5]; // ww

//  printf("uu   = %f\n", uu[0][0]);
//  printf("vv   = %f\n", uu[1][1]);
//  printf("ww   = %f\n", uu[2][2]);

  /* compute eigenvalues. these are not sorted */
  int nrot;
  jacobi(uu, d, vv, DIM, & nrot);

//NEEDED?  /* sort */
//NEEDED?  for(int n=0; n<DIM; n++) {
//NEEDED?    for(int m=0; m<DIM-1; m++) {
//NEEDED?      if(d[m] > d[m+1]) { 
//NEEDED?        real dtmp = d[m]; 
//NEEDED?        d[m] = d[m+1]; 
//NEEDED?        d[m+1] = dtmp;
//NEEDED?      }
//NEEDED?    }
//NEEDED?  }

//DEBUG//printf("d[0] = %f\n", d[0]);
//DEBUG//printf("d[1] = %f\n", d[1]);
//DEBUG//printf("d[2] = %f\n", d[2]);
//DEBUG//printf("nrot = %d\n", nrot);

  /* compute dd (norm?) of d */
  real dd = 0.0;
  for (int i=0; i<DIM; i++) {
    real r = fabs(d[i]);
    d[i] = sqrt(r);
    dd += r;
  }
  dd = sqrt( dd ); 
//DEBUG//printf("dd = %f\n", dd);

  /* compute velocity */
  for (int i=0; i<DIM; i++) v[i]=0.0;

  for (int ie=0; ie<ne; ie++) { 

//DEBUG//printf("ie = %d\n", ie);
    /* compute non-homogeneous k */
    const real SMALL = 1.0e-30;
    for(int i=0; i<DIM; i++) {
      if(d[i] > SMALL) k[i] = K[ie][i] / (tt * d[i]);
      else             k[i] = K[ie][i] / (tt * dd);
//DEBUG//printf("k[i] = %f\n", k[i]);
    }


    /* use non-homogeneous k to compute new velocity */
    a=sclp(k,x)+Omega[ie]*t/tt;
    c=cos(a); s=sin(a);
    for (int i=0; i<DIM; i++) {
      v[i]+=U1[ie][i]*c+U2[ie][i]*s;
//DEBUG//printf("  v[i] = %f\n", v[i]);
    }
  }
  /* anisotropy */
  real fe=sqrt(2./(real)ne);
//DEBUG//printf("fe   = %f\n", fe  );
  for (int i=0; i<DIM; i++) {
    v[i] *= fe*d[i];
//DEBUG//printf("d[i] = %f\n", d[i]);
//DEBUG//printf("v[i] = %f\n", v[i]);
  }

  /* free memory */
  dealloc2d(&uu);
  dealloc2d(&vv);
  dealloc1d(&d);
}

/******************************************************************************/
 main() {
   real x[3],v[3];
   real uu[6];          // Reynolds stresses 
   int nspec=1000;      // spectral sample size
   real turbtime=1.0;   // turbulence time-scale
   real turblength=1.0; // turbulence length-scale
 
   genspec(&nspec);
 
   int  nt=100;   // number of time steps
   real dt=0.5;   // time-step
   int  nx= 32;   // number space-points
   real dx=1.0;   // their separation
   int  ny= 32;   // number space-points
   real dy=1.0;   // their separation
 
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
       for(int iy=0; iy<ny; iy++) 
        {
         real t=dt*it;
         x[0]=dx*ix;
         x[1]=dy*iy;
         x[2]=0.0;
 
         uu[0] = 1.00; // uu
         uu[1] = 1.00; // vv
         uu[2] = 1.00; // ww
         uu[3] = 0.0;  // uv
         uu[4] = 0.0;  // uw
         uu[5] = 0.0;  // vw

         genvec(t,x,turbtime,uu,v);
/////////printf("%f %f %f %f %f\n", t, x[2], v[0], v[1], v[2]);

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
     for(int ix=0; ix<nx; ix++) 
       for(int iy=0; iy<ny; iy++) 
         fprintf(fp, "%f\n", u2d[ix][iy]);
     for(int ix=0; ix<nx; ix++) 
       for(int iy=0; iy<ny; iy++) 
         fprintf(fp, "%f\n", v2d[ix][iy]);
     for(int ix=0; ix<nx; ix++) 
       for(int iy=0; iy<ny; iy++) 
         fprintf(fp, "%f\n", w2d[ix][iy]);
    
     fprintf(fp,"endgmv\n");
    
     fclose(fp);
   }
 }

/*-----------------------------------------------------------------------------+
 '$Id: nhturb_main.cpp,v 1.3 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
