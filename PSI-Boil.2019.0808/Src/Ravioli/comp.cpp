#include "comp.h"

/******************************************************************************/
std::ostream & operator << (std::ostream &ost, const Comp & com) {

  switch(com.val) {
    case( 0): ost << "u "; break;
    case( 1): ost << "v "; break;
    case( 2): ost << "w "; break;
  }

  return ost;
}
