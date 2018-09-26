#include "dispersed.h"

/******************************************************************************/
void Dispersed::check_add_particles() {
/*----------------------------------------------------------------------------+
|  This function is called from the main.cpp for bubble column simulation     |
|  This is to check if the bubbles that enters the domain are injected        |
|  inside bubbles already in the domain (which still are very close from the  |
|  injection points).                                                         |
|  If it turns they will be injected inside already existig bubbles: then     |
|  do not inject them. This doesn't happen so often, but it is crucial to     |
|  check it. Otherwise, if bubble is injected in another bubble, then code    |
|  crash may happen.                                                          |
+----------------------------------------------------------------------------*/

  loop_start_again: 

  cell_link();

  for(int pa=0; pa < size(); pa++) { 

    const real xp = particles[pa].x();
    const real yp = particles[pa].y();
    const real zp = particles[pa].z();
    int ic = int (trunc ( ( xp - dom->global_min_x() ) / diam_x));
    int jc = int (trunc ( ( yp - dom->global_min_y() ) / diam_y));
    int kc = int (trunc ( ( zp - dom->global_min_z() ) / diam_z));

    if(ic == NX_coarse) ic--; 
    if(jc == NY_coarse) jc--;
    if(kc == NZ_coarse) kc--;
    int west   = ic -1; int east  = ic +1;
    int south  = jc -1; int north = jc +1;
    int bottom = kc -1; int top   = kc +1;

    /* take care of boundary conditions */
    if(west   ==       OFF)  west   += 1;
    if(east   == NX_coarse)  east   -= 1;
    if(south  ==       OFF)  south  += 1;
    if(north  == NY_coarse)  north  -= 1;
    if(bottom ==       OFF)  bottom += 1;
    if(top    == NZ_coarse)  top    -= 1;

    for (int ii = west; ii <= east; ii++) {
      for (int jj = south; jj <= north; jj++) {
        for (int kk = bottom; kk <= top; kk++) {

          /* start from the head of the chain */
          int pb = cell[ii][jj][kk];
           
          while (pb > pa) {

            const real dx = particles[pa].x() - particles[pb].x();
            const real dy = particles[pa].y() - particles[pb].y();
            const real dz = particles[pa].z() - particles[pb].z();
            const real dist_ab  = sqrt(dx*dx + dy*dy + dz*dz);
            const real radius_a = particles[pa].d() * 0.5;
            const real radius_b = particles[pb].d() * 0.5;

            real rel_dis = dist_ab - (radius_a + radius_b - boil::pico);
            if(rel_dis < 0.0) {
              /* remove the higher index */
              erase(pb); 
              boil::oout <<"Do Not Add, "   << pb << boil::endl;
              goto loop_start_again; 
            } 

            pb = link[pb];

          } /* while (pb > pa) */
        }
      }
    } 

  } /* for pa */

  delete [] link;
}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed_check_add_particles.cpp,v 1.1 2015/08/19 11:29:51 badreddine Exp $'/ 
+-----------------------------------------------------------------------------*/
