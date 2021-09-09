#include "commonheattransfer.h"
/******************************************************************************/
bool CommonHeatTransfer::edge(const Sign dir, const Comp & m,
                              const int i, const int j, const int k) const {
/***************************************************************************//*** 
*  \brief check if we are at domain edge  
*******************************************************************************/

  if       (m==Comp::i()) {
    if( (i==tpr.si()&&dir<0&&!tpr.bc().type(Dir::imin(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::imin()))
      ||
        (i==tpr.ei()&&dir>0&&!tpr.bc().type(Dir::imax(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::imax())) )
      return true;

  } else if(m==Comp::j()) {
    if( (j==tpr.sj()&&dir<0&&!tpr.bc().type(Dir::jmin(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::jmin()))
      ||
        (j==tpr.ej()&&dir>0&&!tpr.bc().type(Dir::jmax(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::jmax())) )
      return true;
  } else {
    if( (k==tpr.sk()&&dir<0&&!tpr.bc().type(Dir::kmin(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::kmin()))
      ||
        (k==tpr.ek()&&dir>0&&!tpr.bc().type(Dir::kmax(),BndType::periodic())
                           &&!tpr.bc().type_decomp(Dir::kmax())) )
      return true;
  }

  return false;
}
