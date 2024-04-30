#include "simulation_time.h"

/******************************************************************************/
void Times::print_total() const {
/*-----------------------------------------------------------------------------+ 
|  prints information after the creation of the object                         |
+-----------------------------------------------------------------------------*/

  boil::oout << "##############################################" << boil::endl;
  boil::oout << "#                                             " << boil::endl;
  boil::oout << "# TOTAL NUMBER OF TIME STEP: " << total_steps() << boil::endl;
  boil::oout << "#                                             " << boil::endl;
  boil::oout << "# TOTAL SIMULATION TIME:     " << total_time()  << boil::endl;
  boil::oout << "#                                             " << boil::endl;
  boil::oout << "##############################################" << boil::endl;
} 

/******************************************************************************/
void Times::print_current() const {
/*-----------------------------------------------------------------------------+ 
|  prints information about the current time step                              |
+-----------------------------------------------------------------------------*/
  if (print_cinfo) {
    boil::oout << "#############################################" << boil::endl;
    boil::oout << "# STEP: " << current_step() << " / " 
                             << total_steps()                     << boil::endl;
    boil::oout << "# TIME: " << current_time() << " / " 
                             << total_time()                      << boil::endl;
    boil::oout << "# WALL: " << boil::timer.current_min() << " [min]; " 
                             << boil::timer.current_hour() << " [h]"          
                                                                  << boil::endl;
    boil::oout << "#############################################" << boil::endl;
  }
}
