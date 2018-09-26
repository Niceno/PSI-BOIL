#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
int load() {
/*----------------------------------------+
|  loads the data obtained from EasyMesh  |
|  and allocates memory                   |
+----------------------------------------*/

  char   dum[10];
  int    num;
  int    nn, ne, ns, nb;
  FILE * in;
 
  /*------------+
  |             |
  |  node data  |
  |             |
  +------------*/
  name[len+1] = 'n';
  printf("Now loading: %s\n", name);

  if((in=fopen(name, "r"))==NULL) {
   fprintf(stderr, "Cannot load file %s !\n", name);
   return OFF;
  }
 
  fscanf(in, "%d\n", &nn);
  node.resize(nn);

  for(int n=0; n<nn; n++)
    fscanf(in, "%s  %lf %lf  %d\n",
                dum, &node[n].x, &node[n].y, &node[n].mark);
  fclose(in);

  printf("Node data loaded !\n"); fflush(stdout);

  /*---------------+
  |                |
  |  element data  |
  |                |
  +---------------*/
  name[len+1] = 'e';
  printf("Now loading: %s\n", name);

  if((in=fopen(name, "r"))==NULL) {
    fprintf(stderr, "Cannot load file %s !\n", name);
    return OFF;
  }

  fscanf(in, "%d", &ne);
  elem .resize(ne);
  dp_dx.resize(ne);
  dp_dy.resize(ne);
 
  for(int e=0; e<ne; e++) {
    fscanf(in, "%s  %d", dum, &elem[e].n);
    if(elem[e].n==3)
      fscanf(in, "%d %d %d",
                 &elem[e].i,  &elem[e].j,  &elem[e].k);
    if(elem[e].n==4)
       fscanf(in, "%d %d %d %d",
                  &elem[e].i,  &elem[e].j,  &elem[e].k, &elem[e].l);
       fscanf(in, "%lf %lf  %d", &elem[e].x, &elem[e].y, &elem[e].mark);
  }
  fclose(in);
  printf("Element data loaded !\n"); fflush(stdout);

  /*------------+
  |             |
  |  side data  |
  |             |
  +------------*/
  name[len+1] = 's';
  printf("Now loading: %s\n", name);

  if((in=fopen(name, "r"))==NULL) {
    fprintf(stderr, "Cannot load file %s !\n", name);
    return OFF;
  }

  fscanf(in, "%d\n", &ns);
  side.resize(ns);
  SxV .resize(ns);
  SyV .resize(ns);
  SxD .resize(ns);
  SyD .resize(ns);
 
  Flux.resize(ns);

  for(int s=0; s<ns; s++)
    fscanf(in, "%s   %d %d   %d %d   %d\n",
               dum,
               &side[s].c, &side[s].d, 
               &side[s].a, &side[s].b,
               &side[s].mark);
  fclose(in);
  printf("Side data loaded !\n"); fflush(stdout);

  /*---------------------------+
  |                            |
  |  boundary conditions data  |
  |                            |
  +---------------------------*/
  name[len+1] = 'b';
  printf("Now loading: %s\n", name);

  if((in=fopen(name, "r"))==NULL) {
    fprintf(stderr, "Cannot load file %s !\n", name);
    return OFF;
  }

  fscanf(in, "%lf\n", &visc);

  fscanf(in, "%d\n", &nb);
  bound.resize(nb);

  for(int n=0; n<nb; n++) {
    fscanf(in, "%d", &num);
    fscanf(in, "%d  %lf %lf %lf",
                &bound[num].type, 
                &bound[num].u, 
                &bound[num].v, 
                &bound[num].T);
  }
  fclose(in);

  printf("Boundary conditions data loaded !\n"); fflush(stdout);

  /*------------------------------+
  |  allocate the working arrays  |
  +------------------------------*/
  Auvw.con.resize(nn);
  Auvw.val.resize(nn);
  u       .resize(nn);
  v       .resize(nn);
  uo      .resize(nn);
  vo      .resize(nn);
  bu_b    .resize(nn);
  bv_b    .resize(nn);
  bu_Co   .resize(nn);
  bv_Co   .resize(nn);
  bu_Coo  .resize(nn);
  bv_Coo  .resize(nn);
  bu_Do   .resize(nn);
  bv_Do   .resize(nn);
  bu_n    .resize(nn);
  bv_n    .resize(nn);
  bu_P    .resize(nn);
  bv_P    .resize(nn);
  bu      .resize(nn);
  bv      .resize(nn);
 
  Aphi.con.resize(nn);
  Aphi.val.resize(nn);
  phi     .resize(nn);
  bphi    .resize(nn);
  p       .resize(nn);
 
  /*------------------------------+
  |  Count the neighbours of ...  | 
  +------------------------------*/
  for(int s=0; s<ns; s++) {
    node[side[s].c].n++;
    node[side[s].d].n++;
  }

  for(int n=0; n<nn; n++) {
    Auvw.con[n].resize(node[n].n+1);
    Auvw.val[n].resize(node[n].n+1);
  }

  for(int n=0; n<nn; n++) {
    Aphi.con[n].resize(node[n].n+1);
    Aphi.val[n].resize(node[n].n+1);
  }

  return ON;
}

/*-----------------------------------------------------------------------------+
 '$Id: io_load.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/


