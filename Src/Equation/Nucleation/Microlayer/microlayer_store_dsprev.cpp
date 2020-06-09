#include "microlayer.h"
//#define DEBUG

/******************************************************************************/
void Microlayer::store_dSprev() {

#ifdef DEBUG
  std::cout<<"str_dSprev: "<<boil::cart.iam()<<"\n";
#endif

  int kk;
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);
    dSprev[i][j][k] = area_vapor(Sign::neg(),Comp::k(),i,j,k);
    kk=k;
  }
  str_dSprev=true;
}
