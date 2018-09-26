#include "dispersed.h"

/******************************************************************************/
void Dispersed::couple_interface() {
/*-------------------------------------------------+
|  couple dispersed phase with interface tracking  |                                              
+-------------------------------------------------*/

  loop_b_start:

  for_p(p) {
    real vol_dis = 0.0;
    real vol_con = 0.0;
    real r = 0.5 * particles[p].d();
    real delta_x = particles[p].box_dx();
    real delta_y = particles[p].box_dy();
    real delta_z = particles[p].box_dz();

    /* volume of separated field in a particle box, above the surface */
    real v_critical = delta_x * delta_y * (0.5 * delta_z - r);

    real rho_d = flu->rho(dispersed);
    real rho_c = flu->rho(continuous);
 
    for_pijk(p, i, j, k) { /* browse in particle "p" box */

      real sep = (*col)[i][j][k];
      if(sep > 1.0) sep = 1.0;
      if(sep < 0.0) sep = 0.0;
      real d_sep = dispersed_fraction(sep);
      vol_dis += dV(i,j,k) * d_sep;; 
    }

    boil::cart.sum_real(& vol_dis);
    assert(vol_dis >= 0.0);

    /* gas particle in liquid */
    if( rho_d < rho_c ) {
      if(vol_dis > v_critical) {
        merge(p);
        boil::oout << "Merging Gas bubbles in Liquid.... "   
                   << "@dispersed_advance; particle ... " << p  
                   << " entered continuous phase. looping again!" 
                   << boil::endl;
        goto loop_b_start;
      }
    }

    /* liquid (or solid) particle in gas */
    if( rho_d > rho_c ) {
      if(vol_dis > v_critical) {
        merge(p);
        boil::oout << "Merging DROPLETS in Liquid.... "    
                   << "@dispersed_advance; particle... " << p  
                   << " entered continuous phase. looping again!" 
                   << boil::endl;
        goto loop_b_start;
      }
    }
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: dispersed_couple_interface.cpp,v 1.2 2015/08/28 11:32:29 badreddine Exp $'/ 
+-----------------------------------------------------------------------------*/
