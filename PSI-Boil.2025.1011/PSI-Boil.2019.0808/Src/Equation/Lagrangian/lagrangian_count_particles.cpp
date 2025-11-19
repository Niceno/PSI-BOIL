#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::count_particles() {
/*-----------------------------------------------------+
|  telling how many particles in a simulation we have  |
+-----------------------------------------------------*/

  /* the disctinciont below is only important for the messages being printed */

  assert(lagrangian == 0.56789);  //mark

  // particle
  if(size() > 0) {
    boil::oout << "@lagrangian_advance; advancing " 
               << size() << " particle." << boil::endl;
  } else {
    boil::oout << "@lagrangian_advance; no particle to advance. " 
               << boil::endl;
    boil::oout << "-------------end-------------" << boil::endl;
    getchar();
    getchar();
    getchar();
  }
  return;

}

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_count_particles.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
