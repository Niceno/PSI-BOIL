#include "floodfill.h"

/***************************************************************************//**
*  Compares color function values on either side of decomposed domain.
*******************************************************************************/
void Floodfill::compare_rgn(const int & rgnin, const int & rgnout) {
  if (rgnin * rgnout > 0) {  
    int i=0;
    for (i; i<rgnid_match_arrayi.size(); i++) {
      if ((rgnid_match_arrayi[i] == rgnin) && 
          (rgnid_match_arrayj[i] == rgnout))  {
        break; 
      }
    }
    if (i == rgnid_match_arrayi.size()) {
      rgnid_match_arrayi.push_back(rgnin); 
      rgnid_match_arrayj.push_back(rgnout); 
    }
  }
}
