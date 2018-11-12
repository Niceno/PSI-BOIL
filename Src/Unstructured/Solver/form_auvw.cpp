#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "definitions.h"
#include "variables.h"

/******************************************************************************/
void form_Auvw() {

  /*----------------+
  |  viscous terms  |
  +----------------*/
  for_sides(s) {

    const int c=side[s].c;
    const int d=side[s].d;

    /* factor 0.5 in front of diffusive terms is for semi-implicit treatment */
    Auvw.val[c][side[s]._d] = - 0.5 * visc * side[s].l / side[s].h;  // 0.5 <- v1
    Auvw.val[d][side[s]._c] = - 0.5 * visc * side[s].l / side[s].h;  // 0.5 <- v1

    Auvw.val[c][0] -= Auvw.val[c][side[s]._d];
    Auvw.val[d][0] -= Auvw.val[d][side[s]._c];
  }

  /*-----------------+
  |  unsteady terms  |
  +-----------------*/
  for_nodes(n)
    if(node[n].mark==0)
      Auvw.val[n][0] += node[n].Area/dt;

  /*-----------------------------------+
  |  take care of boundary conditions  |
  +-----------------------------------*/
  for_sides(s) {
    const int c=side[s].c;
    const int d=side[s].d;

/*-----Inlet. or adiabatic wall-----*/    // more sophisticated b.c. checking v3
    if(bound[node[c].mark].type==GEO || bound[node[c].mark].type==NAT  ||
       bound[node[d].mark].type==GEO || bound[node[d].mark].type==NAT) {

      if(node[c].mark!=0) {
        Auvw.val[c][0]          = 1.0;
        Auvw.val[c][side[s]._d] = 0.0;
        Auvw.val[d][side[s]._c] = 0.0;
        bu_b[c] = bound[node[c].mark].u; 
        bv_b[c] = bound[node[c].mark].v;
      } else {
        bu_b[c] += 0.5 * visc * side[s].l / side[s].h * bound[node[d].mark].u;  
        bv_b[c] += 0.5 * visc * side[s].l / side[s].h * bound[node[d].mark].v;
      } 

      if(node[d].mark!=0 ) {
        Auvw.val[d][0]          = 1.0;
        Auvw.val[d][side[s]._c] = 0.0;
        Auvw.val[c][side[s]._d] = 0.0;
        bu_b[d] = bound[node[d].mark].u;
        bv_b[d] = bound[node[d].mark].v; 
      } else {
        bu_b[d] += 0.5 * visc * side[s].l / side[s].h * bound[node[c].mark].u;  
        bv_b[d] += 0.5 * visc * side[s].l / side[s].h * bound[node[c].mark].v;
      } 

    }
  } /* through sides */
}
