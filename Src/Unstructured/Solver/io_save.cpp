#include <cstdio>
#include <cstdlib> 
#include <cstring>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
int save() {
/*-------------------------+
|  create gmv result file  |
+-------------------------*/

  FILE * out;

  name[len+1] = 'r';

  name[len+2] = '\0';
  strcat(name, ".gmv");
  if((out=fopen(name, "w"))==NULL) {
    fprintf(stderr, "Cannot save file !\n");
    return OFF;
  }

  /*--------+
  |  start  | 
  +--------*/
  fprintf(out, "gmvinput ascii\n");

  /*--------+
  |  nodes  | 
  +--------*/
  fprintf(out, "nodes %d\n", (int)node.size());

  for_nodes(n) 
    fprintf(out, "%lf\n", node[n].x);

  for_nodes(n) 
    fprintf(out, "%lf\n", node[n].y);

  for_nodes(n) 
    fprintf(out, "0.0\n");
 
  /*--------+
  |  cells  | 
  +--------*/
  fprintf(out, "cells %d\n", (int)elem.size());
  for_elems(e) { 
    if(elem[e].n==3)
      fprintf(out, "tri 3 %d %d %d\n", 
        elem[e].i+1, elem[e].j+1, elem[e].k+1);

    if(elem[e].n==4)
      fprintf(out, "quad 4 %d %d %d %d\n", 
        elem[e].i+1, elem[e].j+1, elem[e].k+1, elem[e].l+1);
  }

  /*-----------+
  |  solution  | 
  +-----------*/
  fprintf(out, "velocity 1\n");

  for_nodes(n) 
    fprintf(out, "%lf\n", u[n]);

  for_nodes(n) 
    fprintf(out, "%lf\n", v[n]);

  for_nodes(n) 
    fprintf(out, "%lf\n", 0.0);

  fprintf(out, "variable\n");

  fprintf(out, "p 1\n");
  for_nodes(n) 
    fprintf(out, "%9.6e\n", p[n]);

   fprintf(out, "endvars\n");

  /*------+
  |  end  | 
  +------*/
  fprintf(out, "endgmv\n");

  fclose(out);

  return ON;
}

/*-----------------------------------------------------------------------------+
 '$Id: io_save.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
