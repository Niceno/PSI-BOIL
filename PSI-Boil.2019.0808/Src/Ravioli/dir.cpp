#include "dir.h"

/***************************************************************************//**
*  Prints the name of the current limiter. Needed only for information or 
*  debugging, if ever at all. 
*******************************************************************************/
std::ostream & operator << (std::ostream &ost, const Dir & dir) {

  switch(dir.val) {
    case(-1): ost << "undefined "; break;
    case( 0): ost << "imin ";      break;
    case( 1): ost << "imax ";      break;
    case( 2): ost << "jmin ";      break;
    case( 3): ost << "jmax ";      break;
    case( 4): ost << "kmin ";      break;
    case( 5): ost << "kmax ";      break;
    case( 6): ost << "ibody ";     break;
  }
  
  return ost;
}
