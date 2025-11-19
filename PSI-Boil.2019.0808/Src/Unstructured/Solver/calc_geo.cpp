#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
void calc_geo() {

/*-----------+
|            |
|  Find ...  |                       
|            |
+-----------*/
 for_sides(s) {
   const real xc=node[side[s].c].x; 
   const real yc=node[side[s].c].y;
   const real xd=node[side[s].d].x; 
   const real yd=node[side[s].d].y;

/*-------------------------+
|  ... triangle sides (h)  |
+-------------------------*/
   side[s].h = sqrt( (xc-xd)*(xc-xd) + (yc-yd)*(yc-yd) );

/*--------------------------------+
|  ... Voronoi polygon sides (l)  | <== ovdje je mozda bila greska !!
+--------------------------------*/
   real xa, ya, xb, yb;

   if(side[s].a != OFF)
    {xa=elem[side[s].a].x; ya=elem[side[s].a].y;}
   else
    {xa=0.5*(xc+xd); ya=0.5*(yc+yd);}

   if(side[s].b != OFF)
    {xb=elem[side[s].b].x; yb=elem[side[s].b].y;}
   else
    {xb=0.5*(xc+xd); yb=0.5*(yc+yd);}

   side[s].l = sqrt( (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) );
  
/*--------------------------------------------+
|  Voronoi nodes ... Delaunay triangle areas  |
+--------------------------------------------*/
   if(side[s].a != OFF)
     elem[side[s].a].Area+=.5*((xd-xa)*(yd+ya)+(xa-xc)*(ya+yc)+(xc-xd)*(yc+yd));
   if(side[s].b != OFF)
     elem[side[s].b].Area-=.5*((xd-xb)*(yd+yb)+(xb-xc)*(yb+yc)+(xc-xd)*(yc+yd));

/*------------------------------------------+
|  ... components of Voronoi polygon sides  |
|  It is important to note that the sides   |
|       are oriented from c to d node       |
+------------------------------------------*/
   SxV[s] = ya - yb;
   SyV[s] = xb - xa;
   SxD[s] = yc - yd;
   SyD[s] = xd - xc;

/*----------------------------+
|  ... Voronoi polygon areas  |
+----------------------------*/
   node[side[s].c].Area += 0.25 * side[s].h * side[s].l;
   node[side[s].d].Area += 0.25 * side[s].h * side[s].l;
  }

/*%%%%%%%%%
%  Check  %
%%%%%%%%%*/
 {
  real Ar;

  Ar=0.0;
  for_nodes(n) Ar+=node[n].Area;
  printf("Area  I = %lf\n", Ar);

  Ar=0.0;
  for_elems(e) Ar+=elem[e].Area;
  printf("Area II = %lf\n", Ar);

 }

/*---------------------------------+
|  find the interpolation factors  |            // new block from v4
|      for pressure gradients      |
+---------------------------------*/
 {
 /*----- find the number of elements surrounding the node -----*/
  for_elems(e) {
    node[elem[e].i].Ne++;
    node[elem[e].j].Ne++;
    node[elem[e].k].Ne++;
    if(elem[e].n==4) 
      node[elem[e].l].Ne++;
  }

   /*----- allocate memory for storing the interpolation factrs -----*/
  for_nodes(n) {
    node[n].elem.resize( node[n].Ne );
    node[n].fe  .resize( node[n].Ne );
    node[n].Ne   = 0;  
  }

  for_elems(e) {
 
    /* i */
    const int i=elem[e].i;
    node[i].Ne++;
    node[i].elem[node[i].Ne-1] = e;
    node[i].fe[node[i].Ne-1]   = 1.0/elem[e].Area;

    /* j */
    const int j=elem[e].j;
    node[j].Ne++;
    node[j].elem[node[j].Ne-1] = e;
    node[j].fe[node[j].Ne-1]   = 1.0/elem[e].Area;

    /* k */
    const int k=elem[e].k;
    node[k].Ne++;
    node[k].elem[node[k].Ne-1] = e;
    node[k].fe[node[k].Ne-1]   = 1.0/elem[e].Area;

    /* l */
    if(elem[e].n==4)
     {const int l=elem[e].l;
      node[l].Ne++;
      node[l].elem[node[l].Ne-1] = e;
      node[l].fe[node[l].Ne-1]   = 1.0/elem[e].Area;}
   }

  /* normalize the interpolation factors */
  for_nodes(n) {
   
    real sum=0;
    for(int e=0; e<node[n].Ne; e++)
      sum += node[n].fe[e];

    for(int e=0; e<node[n].Ne; e++)
      node[n].fe[e] /= sum;
  }

  /*%%%%%%%%%
  %  Check  %
  %%%%%%%%%*/
  for_nodes(n) {

    printf("%3d: %3d ", n, node[n].Ne);  
    for(int e=0; e<node[n].Ne; e++)
      printf("%lf  ", node[n].fe[e]);
    printf("\n");
   } 

 } /* end of block */

}
