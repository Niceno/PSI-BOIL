/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%  Notes:                                                                      %
%                                                                              %
%  New in this version are the interpolation factors for pressure gradients    %
%  computed for each node.  That required a more elaborate Node structure.     %
%                                                                              %
%  Logistic detail: load and save functions have been moved to "utils.h"       %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include "structures.h"
#include "definitions.h"

/* global variable definitions */
std::vector<struct ele> elem;
std::vector<struct sid> side;
std::vector<struct nod> node;
std::vector<struct bc>  bound;

real dt, visc, CFL;  /* what is this for ? */ 
char   name[256]; 
int    len;
std::vector<real> SxV, SyV, SxD, SyD;

class matrix Auvw, Aphi;
std::vector<real> u, v, uo, vo, 
                  bu_b,   bv_b,   /* boundary conditions */
                  bu_Co,  bv_Co,  /* convective terms */
                  bu_Coo, bv_Coo, /* convective terms */
                  bu_Do,  bv_Do,  /* diffusuvwe terms */
                  bu_n,   bv_n,   /* nonstationary terms */
                  bu_P,   bv_P,   /* pressure terms */
                  bu,     bv,     /* total fluxes */
                  phi, bphi, p, Flux, 
                  dp_dx, dp_dy;

/* global function definitions */
int load();
int save();
int solve_cg(const matrix            & A, 
                   std::vector<real> & x, 
             const std::vector<real> & b, 
             const real                tol, 
             const int                 prec, 
             const int                 niter);
int solve_bicg(const matrix            & A, 
                     std::vector<real> & x, 
               const std::vector<real> & b, 
               const real                tol, 
               const int                 prec, 
               const int                 niter);
void calc_geo();
void calc_con();
void form_Auvw();
void form_buvw();
void form_Aphi();
int  form_bphi();
int  project_uvw();

/******************************************************************************/
main(int argc, char *argv[]) {

  printf("%d\n", argc);

  if(argc==2) {
    strcpy(name, argv[1]);
    len=strlen(name);
    if(name[len-2]=='.')
      if(name[len-1]=='d' || name[len-1]=='D' )
        name[len-2]='\0';
  }
  strcat(name, ".x");
  printf("ime: %s\n", name);

  /*----------------+
  |  load the mesh  |
  +----------------*/
  if(load()!=ON)
    return OFF;

  /*-----------------------------------+
  |  calculate geometrical quantities  |
  +-----------------------------------*/
  calc_geo();
  printf("calc geo !\n"); fflush(stdout);
  calc_con();

  printf("Time step ?");
  scanf("%lf", &dt);
  form_Aphi();
  form_Auvw();

  int Ndt;
  do {
    printf("Number of time steps ?");
    scanf("%d", &Ndt);

    for(int n=0; n<Ndt; n++) {
      printf("%4d: ", n);

      form_buvw();
      printf(" u:%d ", solve_bicg(Auvw, u, bu, 1e-10, 2, node.size()));
      printf(" v:%d ", solve_bicg(Auvw, v, bv, 1e-10, 2, node.size()));
   
      form_bphi();
   
      printf(" phi:%d",     solve_cg(Aphi, phi, bphi, 1e-10, 2, node.size()));
      printf(" CFL:%5.3lf", CFL);
 
      project_uvw();

      printf("\n");
    }

   save();
  } while(Ndt > 0);

  return ON;
}

/*-----------------------------------------------------------------------------+
 '$Id: main.cpp,v 1.2 2015/09/01 09:51:02 niceno Exp $'/
+-----------------------------------------------------------------------------*/
