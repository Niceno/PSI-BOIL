#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::advance(Vector * force) {
/*==========================================================================+
|                                                                           |
|              the main algorithm to advance the lagrangian phase            |
|                                                                           |
|===========================================================================|
|                                                                           |
|  1) check_particles(): check the number of particles that exist in the    |
|     domain.                                                               |
|                                                                           |
|  2) collisions(): handle partilce-particle and particle-wall collisions.  |
|                                                                           |
|  3) box_particles(): create a box around every particle.                  | 
|                                                                           |
|  4) couple_interface(): check if the particles are close to the free      |
|     surface and if yes merge them.                                        |
|                                                                           |
|  5) repaint(): update the particles'color function.                       | 
|                                                                           |
|  6) forces(): compute the forces acting on particles: buoyancy, drag,     | 
|     added mass, lift, wall lubrication force, and then update             |
|     particles's velocity.                                                 |           
|                                                                           |
|  7) volume_averaged(): calculate volume-averaged velocity and density     |
|     over particle's volume.                                               |        
|                                                                           |
|  8) modeled(force): calculate the modeled_force that will set the         | 
|     volume-averaged velocity to the particle's calculated velocity.       | 
|                                                                           |
+==========================================================================*/

  boil::timer.start("lagrangian advance");

  count_particles();

  collisions();  //update-particle-position-xyz

  box_particles();

  //couple_interface();  //touch-then-end

  /*#if CIT == false
    repaint(); 
  #endif*/

  forces();  //update-particle-velocity-uvw

  //volume_averaged();  //preparations-for-interphase-coupling

  //modeled(force);  //two-way-coupling

  boil::timer.stop("lagrangian advance");

  /*OPR("lagrangian_advance.cpp");
  getchar();*/

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_advance.cpp,v 1.16 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
