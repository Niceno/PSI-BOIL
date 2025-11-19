#include "dispersed.h"

/******************************************************************************/
void Dispersed::count_particles() {
/*-----------------------------------------------------+
|  telling how many particles in a simulation we have  |
+-----------------------------------------------------*/

  /* the disctinciont below is only important for the messages being printed */

  /* bubble */
  if( dispersed < continuous ) {  
    if(size() > 0) {
      boil::oout << "@dispersed_advance; advancing " 
                 << size() << " bubbles." << boil::endl;
    } else {
      boil::oout << "@dispersed_advance; no bubbles to advance. " 
                 << boil::endl;
    }
    return;

  /* droplet */
  } else { 
    if(size() > 0) {
      boil::oout << "@dispersed_advance; advancing " 
                 << size() << " droplets." << boil::endl;
    } else {
      boil::oout << "@dispersed_advance; no droplets to advance. " 
                 << boil::endl;
    }
     return;
  }

}
