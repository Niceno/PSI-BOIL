#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "variables.h"
#include "definitions.h"

/******************************************************************************/
int form_bphi() {

  /*------------------------------+
  |  erase the source term first  |
  +------------------------------*/
  for_nodes(n) {
    phi[n] =0.0;
    bphi[n]=0.0;
  }

  /*-----------------------+
  |  fill the source term  |
  +-----------------------*/
  CFL=0.0;

  /*-----Fluxes inside the domain-----*/ 
  for_sides(s) {
    const int c = side[s].c;
    const int d = side[s].d;

    Flux[s] = 0.5 * (  (  u[c]-bu_P[c]/Auvw.val[c][0]
                         +u[d]-bu_P[d]/Auvw.val[d][0] ) * SxV[s] 
                                      + 
                       (  v[c]-bv_P[c]/Auvw.val[c][0]
                         +v[d]-bv_P[d]/Auvw.val[d][0] ) * SyV[s] );
  
    Flux[s] += (p[c]-p[d])/side[s].h * side[s].l * dt;
    /*         |                   |   |       |
               |<----- dp/dn ----->|   |<--l-->|
     */  
    CFL = boil::maxr(CFL, fabs(Flux[s]/side[s].l*dt/side[s].h));

    bphi[d] += Flux[s];
    bphi[c] -= Flux[s];
  }

  /*-----Boundary fluxes-----*/
  real InFlow  = 0.0;
  real OutFlow = 0.0;

  for_sides(s) {
    const int a = side[s].a;
    const int b = side[s].b;
    const int c = side[s].c;
    const int d = side[s].d;

    if(side[s].mark != 0) {
      real Flux_c = 0.0;
      real Flux_d = 0.0;
  
      const real xs = 0.5*(node[c].x + node[d].x);
      const real ys = 0.5*(node[c].y + node[d].y);

      if(side[s].a==OFF) {
        if(ys < elem[b].y) 
         {Flux_c += 0.25*fabs(SyD[s]) * (v[c]+v[d]);            // Flux_c from v3
          Flux_d += 0.25*fabs(SyD[s]) * (v[c]+v[d]);}           // Flux_d from v3
        else      
         {Flux_c -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);            // Flux_c from v3
          Flux_d -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);}           // Flux_d from v3
        if(xs < elem[b].x) 
         {Flux_c += 0.25*fabs(SxD[s]) * (u[c]+u[d]);            // Flux_c from v3
          Flux_d += 0.25*fabs(SxD[s]) * (u[c]+u[d]);}           // Flux_d from v3
        else      
         {Flux_c -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);            // Flux_c from v3
          Flux_d -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);}           // Flux_d from v3
      }

      if(side[s].b==OFF) {
        if(ys < elem[a].y) 
         {Flux_c += 0.25*fabs(SyD[s]) * (v[c]+v[d]);            // Flux_c from v3
          Flux_d += 0.25*fabs(SyD[s]) * (v[c]+v[d]);}           // Flux_d from v3
        else      
         {Flux_c -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);            // Flux_c from v3
          Flux_d -= 0.25*fabs(SyD[s]) * (v[c]+v[d]);}           // Flux_d from v3
        if(xs < elem[a].x) 
         {Flux_c += 0.25*fabs(SxD[s]) * (u[c]+u[d]);            // Flux_c from v3
          Flux_d += 0.25*fabs(SxD[s]) * (u[c]+u[d]);}           // Flux_d from v3
        else      
         {Flux_c -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);            // Flux_c from v3
          Flux_d -= 0.25*fabs(SxD[s]) * (u[c]+u[d]);}           // Flux_d from v3
      }

      bphi[c]+=Flux_c;                                                      // v3
      bphi[d]+=Flux_d;                                                      // v3

      if( bound[side[s].mark].type == OUT )                                 // v3
        OutFlow += (Flux_c+Flux_d);                                         // v3
      else                                                                  // v3
        InFlow  += (Flux_c+Flux_d);                                         // v3
    }
  }

  /*----------------+
  |  compute error  |
  +----------------*/
  real error = 0.0;
  for_nodes(n) {

    error = boil::maxr(error, fabs(bphi[n]));
    bphi[n] /= dt;
  }
  printf(" m: %8.4e ", error);
 
  return ON;
}

/*-----------------------------------------------------------------------------+
 '$Id: form_bphi.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
