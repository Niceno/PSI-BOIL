#include <cmath>
#include <cstdio>
#include <cstdlib> 

#include "definitions.h"
#include "variables.h"

/******************************************************************************/
void form_buvw() {

  /*-----------------+
  |  unsteady terms  |
  +-----------------*/
  for_nodes(n)
    if(node[n].mark==0) {
      uo[n] = u[n];
      vo[n] = v[n];
      bu_n[n] = node[n].Area * uo[n]/dt;
      bv_n[n] = node[n].Area * vo[n]/dt;
    }

  /*------------------+
  |  advective terms  | -> boundary not treated yet
  +------------------*/
  for_nodes(n)
    if(node[n].mark==0) {
      bu_Coo[n] = bu_Co[n];
      bv_Coo[n] = bv_Co[n];
      bu_Co[n] = 0.0;
      bv_Co[n] = 0.0;
    }

 for_sides(s) {
   const int c = side[s].c;
   const int d = side[s].d;

   const real u_ = 0.5 * (u[c]+u[d]);
   const real v_ = 0.5 * (v[c]+v[d]);
  
   bu_Co[d] += Flux[s]*u_;
   bu_Co[c] -= Flux[s]*u_;
   bv_Co[d] += Flux[s]*v_;
   bv_Co[c] -= Flux[s]*v_;
 }

  /*----------------+
  |  viscous terms  |
  +----------------*/
  for_nodes(n) {
    bu_Do[n] = 0.0;
    bv_Do[n] = 0.0;
  }

  for_sides(s) {
    const int c=side[s].c;
    const int d=side[s].d;
  
    bu_Do[c] += 0.5 * visc * (u[d]-u[c]) * side[s].l / side[s].h;
    bu_Do[d] += 0.5 * visc * (u[c]-u[d]) * side[s].l / side[s].h;
    bv_Do[c] += 0.5 * visc * (v[d]-v[c]) * side[s].l / side[s].h;
    bv_Do[d] += 0.5 * visc * (v[c]-v[d]) * side[s].l / side[s].h;
  }

  /*-----------------+
  |  pressure terms  | 
  +-----------------*/
  
  /* initialize the pressure terms */
  for_nodes(n) {
    bu_P[n]=0.0;
    bv_P[n]=0.0;
  }

  /* find the pressure derivations on dual complex
     and average them over nodes */
  for_elems(e) {

    if(elem[e].n==3) {
      const int i = elem[e].i; 
      const int j = elem[e].j; 
      const int k = elem[e].k;

      dp_dx[e] = - 0.5 * (  (p[j]+p[i])*(node[j].y-node[i].y) +
                            (p[k]+p[j])*(node[k].y-node[j].y) + 
                            (p[i]+p[k])*(node[i].y-node[k].y) ) / elem[e].Area;  

      dp_dy[e] =   0.5 * (  (p[j]+p[i])*(node[j].x-node[i].x) +
                            (p[k]+p[j])*(node[k].x-node[j].x) + 
                            (p[i]+p[k])*(node[i].x-node[k].x) ) / elem[e].Area;  
    }
    if(elem[e].n==4) {
      const int i = elem[e].i; 
      const int j = elem[e].j; 
      const int k = elem[e].k; 
      const int l = elem[e].l;

      dp_dx[e] = - 0.5 * (  (p[j]+p[i])*(node[j].y-node[i].y) +
                            (p[k]+p[j])*(node[k].y-node[j].y) + 
                            (p[l]+p[k])*(node[l].y-node[k].y) + 
                            (p[i]+p[l])*(node[i].y-node[l].y) ) / elem[e].Area;  

      dp_dy[e] =   0.5 * (  (p[j]+p[i])*(node[j].x-node[i].x) +
                            (p[k]+p[j])*(node[k].x-node[j].x) + 
                            (p[l]+p[k])*(node[l].x-node[k].x) + 
                            (p[i]+p[l])*(node[i].x-node[l].x) ) / elem[e].Area;  
    }
  } /* e */ 

/*----- integrate pressure derivatinos over cells -----*/
  for_nodes(n) {
    for(int e=0; e<node[n].Ne; e++) {
      bu_P[n] += dp_dx[node[n].elem[e]] * node[n].fe[e];
      bv_P[n] += dp_dy[node[n].elem[e]] * node[n].fe[e];
    }
  } 

  for_nodes(n) {
    bu_P[n] *= node[n].Area;
    bv_P[n] *= node[n].Area;
  }

  /*------------------------+
  |  summ of all the terms  |
  +------------------------*/
  for_nodes(n)
    if(node[n].mark==0) {
      bu[n] = bu_b[n] 
            + bu_n[n] 
            + bu_P[n] 
            + bu_Do[n] 
            + (1.5*bu_Co[n] - 0.5*bu_Coo[n]); 
      bv[n] = bv_b[n] 
            + bv_n[n] 
            + bv_P[n] 
            + bv_Do[n] 
            + (1.5*bv_Co[n] - 0.5*bv_Coo[n])
            + 0.5 * (2.0 - node[n].x) * node[n].Area;
    }
}

/*-----------------------------------------------------------------------------+
 '$Id: form_buvw.cpp,v 1.1 2015/09/01 09:41:25 niceno Exp $'/
+-----------------------------------------------------------------------------*/
