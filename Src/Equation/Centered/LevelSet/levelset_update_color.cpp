#include "levelset.h"

/******************************************************************************/
void LevelSet::update_color(Scalar & sca) {
/***************************************************************************//**
*  \brief update color function
*         input & output: sca
*******************************************************************************/

  for_aijk(i,j,k){
    if(phi[i][j][k]<-1.5*dxmin){
      sca[i][j][k]=0.0;
    } else if(phi[i][j][k]>1.5*dxmin){
      sca[i][j][k]=1.0;
    } else {
      sca[i][j][k]=0.5+phi[i][j][k]/(3.0*dxmin)
                  +1.0/(2.0*pi)*sin(pi*phi[i][j][k]/(1.5*dxmin));
    }
  }
  insert_bc_color(sca);
  sca.exchange_all();
  return;
}
