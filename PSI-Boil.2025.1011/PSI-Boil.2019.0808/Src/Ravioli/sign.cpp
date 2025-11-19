#include "sign.h"

/***************************************************************************//**
*  Prints the name of the current sign. Needed only for information or 
*  debugging, if ever at all. 
*******************************************************************************/
std::ostream & operator << (std::ostream &ost, const Sign & sig) {

  switch(sig.val) {
    case( 0): ost << "undefined "; break;
    case(-1): ost << "negative "; break;
    case(+1): ost << "positive "; break;
  }
  
  return ost;
}
