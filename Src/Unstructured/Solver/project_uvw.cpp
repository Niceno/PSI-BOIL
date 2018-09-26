#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
int project_uvw() {

  /*---------------------+
  |  correct the fluxes  |
  |  and the velocities  |
  +---------------------*/
  for_sides(s) {

   const int c = side[s].c;
   const int d = side[s].d;

   const real Flux_ = (phi[c]-phi[d]) * side[s].l / side[s].h*dt;

   const real u_ = (Flux_/side[s].l) * (node[d].x-node[c].x) / side[s].h;
   const real v_ = (Flux_/side[s].l) * (node[d].y-node[c].y) / side[s].h;
              /*   |               |   |                                |
                   |<-----uv'----->|   |<---------- cos & sin --------->|
              */

   Flux[s] += (phi[c]-phi[d])*side[s].l/side[s].h*dt;

   if(node[c].mark==0) u[c] += u_/node[c].n;
   if(node[c].mark==0) v[c] += v_/node[c].n;
   if(node[d].mark==0) u[d] += u_/node[d].n;
   if(node[d].mark==0) v[d] += v_/node[d].n;
  }


  /*----------------------------+
  |  update the pressure field  |
  +----------------------------*/
  for_nodes(n)
    p[n] += phi[n]*dt;

  /*------------------------------+
  |  erase the source term first  |
  +------------------------------*/
   for_nodes(n)
     bphi[n] = 0.0;

  /*-----------------------+
  |  fill the source term  |
  +-----------------------*/
  for_sides(s) {
    const int c = side[s].c;
    const int d = side[s].d;
  
    bphi[d] += Flux[s];
    bphi[c] -= Flux[s];
  }

  for_sides(s) {
    const int a = side[s].a;
    const int b = side[s].b;
    const int c = side[s].c;
    const int d = side[s].d;

    if(side[s].mark!=0) {
      const real xs = 0.5*(node[c].x+node[d].x);
      const real ys = 0.5*(node[c].y+node[d].y);

      if(side[s].a==OFF) {
        if(ys<elem[b].y) 
         {bphi[c] += 0.25*fabs(SyD[s]) * (v[c]+v[d]);
          bphi[d] += 0.25*fabs(SyD[s]) * (v[c]+v[d]);}
        else      
         {bphi[c] -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);
          bphi[d] -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);}
        if(xs<elem[b].x) 
         {bphi[c] += 0.25*fabs(SxD[s]) * (u[c]+u[d]);
          bphi[d] += 0.25*fabs(SxD[s]) * (u[c]+u[d]);}
         else      
         {bphi[c] -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);
          bphi[d] -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);}
      }
      if(side[s].b==OFF) {
        if(ys<elem[a].y) 
         {bphi[c] += 0.25*fabs(SyD[s]) * (v[c]+v[d]);
          bphi[d] += 0.25*fabs(SyD[s]) * (v[c]+v[d]);}
        else      
         {bphi[c] -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);
          bphi[d] -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);}
        if(xs<elem[a].x) 
         {bphi[c] += 0.25*fabs(SxD[s]) * (u[c]+u[d]);
          bphi[d] += 0.25*fabs(SxD[s]) * (u[c]+u[d]);}
        else      
         {bphi[c] -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);
          bphi[d] -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);}
      }
    }
  }

  real error=0.0;
  for_nodes(n)
    error = boil::maxr(error, fabs(bphi[n]));
  printf(" m: %8.4e ", error);

  return ON;
}

/*-----------------------------------------------------------------------------+
 '$Id: project_uvw.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/

