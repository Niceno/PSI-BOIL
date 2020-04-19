#include "phasechange4.h"
/******************************************************************************/
bool PhaseChange4::edge(const Sign dir, const Comp & m,
                        const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief check if we are at domain edge  
*******************************************************************************/

  if       (m==Comp::i()) {
    if( (i==si()&&dir<0&&!tpr.bc().type(Dir::imin(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::imin()))
      ||
        (i==ei()&&dir>0&&!tpr.bc().type(Dir::imax(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::imax())) )
      return true;

  } else if(m==Comp::j()) {
    if( (j==sj()&&dir<0&&!tpr.bc().type(Dir::jmin(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::jmin()))
      ||
        (j==ej()&&dir>0&&!tpr.bc().type(Dir::jmax(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::jmax())) )
      return true;
  } else {
    if( (k==sk()&&dir<0&&!tpr.bc().type(Dir::kmin(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::kmin()))
      ||
        (k==ek()&&dir>0&&!tpr.bc().type(Dir::kmax(),BndType::periodic())
                       &&!tpr.bc().type_decomp(Dir::kmax())) )
      return true;
  }

  return false;
}
