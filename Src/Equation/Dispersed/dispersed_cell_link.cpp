#include "dispersed.h"

/******************************************************************************/
void Dispersed::cell_link() {

  boil::timer.start("dispersed cell_link");

  /* declare the link vector */
  link = new int [size()];  

  /* initialize cell matrix to OFF */
  for(int ic = 0; ic < NX_coarse; ic++) 
    for(int jc = 0; jc < NY_coarse; jc++) 
      for(int kc = 0; kc < NZ_coarse; kc++) 
        cell[ic][jc][kc] = OFF;
   
  /*-------------------------------------------------------------------------+
  |  Browse through all the particles, and fill link and cell. At the end,   |
  |  cell[ic][jc][kc] will contain the higher index of all the particles     |
  |  inside this cell. If there is no particle in this cell, then its value  |
  |  will be OFF == -1.                                                      |     
  +-------------------------------------------------------------------------*/

  for_p(p) {
    const real xp = particles[p].x();
    const real yp = particles[p].y();
    const real zp = particles[p].z();
    int ic = int ( trunc (( xp - dom->global_min_x() ) / diam_x));
    int jc = int ( trunc (( yp - dom->global_min_y() ) / diam_y));
    int kc = int ( trunc (( zp - dom->global_min_z() ) / diam_z));

    if(ic == NX_coarse) ic--;
    if(jc == NY_coarse) jc--; 
    if(kc == NZ_coarse) kc--;

    link[p] = cell[ic][jc][kc];
    cell[ic][jc][kc] = p;

  }   

  boil::timer.stop("dispersed cell_link");
}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed_cell_link.cpp,v 1.1 2015/08/19 11:29:51 badreddine Exp $'/
+-----------------------------------------------------------------------------*/
