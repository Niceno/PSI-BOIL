#include "vof.h"

/***************************************************************************//**
*  \brief Set flag for indicating near-interface region; used in enthalpy and
*         phase change, passed through topology.
*  iflag = -3 : color function <  clrsurf & no interface nearby
*  iflag = +3 : color function >= clrsurf & no interface nearby
*  iflag = -2 : color function <  clrsurf & interface nearby
*  iflag = +2 : color function >= clrsurf & interface nearby
*  iflag = -1 : color function <  clrsurf & at interface
*  iflag = +1 : color function >= clrsurf & at interface 
*  iflag = -1000 solid cell
*******************************************************************************/
void VOF::interfacial_flagging(Scalar & scp) {

  for_aijk(i,j,k) {
    iflag[i][j][k]=-3;
    if(scp[i][j][k]>=phisurf)
      iflag[i][j][k]=3;
  }

  /* i-direction */
  for(int i=si()-1; i<=ei(); i++){
    for_jk(j,k){
      if(Interface(Sign::pos(),Comp::i(),i,j,k)) {
        if(scp[i][j][k]<phisurf) {
          iflag[i  ][j][k]=-2;
          iflag[i+1][j][k]=+2;
        } else {
          iflag[i  ][j][k]=+2;
          iflag[i+1][j][k]=-2;
        }
      }
    }
  }
  /* j-direction */
  for(int j=sj()-1; j<=ej(); j++){
    for_ik(i,k){
      if(Interface(Sign::pos(),Comp::j(),i,j,k)) {
        if(scp[i][j][k]<phisurf) {
          iflag[i][j  ][k]=-2;
          iflag[i][j+1][k]=+2;
        } else {
          iflag[i][j  ][k]=+2;
          iflag[i][j+1][k]=-2;
        }
      }
    }
  }
  /* k-direction */
  for(int k=sk()-1; k<=ek(); k++){
    for_ij(i,j){
      if(Interface(Sign::pos(),Comp::k(),i,j,k)) {
        if(scp[i][j][k]<phisurf) {
          iflag[i][j][k  ]=-2;
          iflag[i][j][k+1]=+2;
        } else {
          iflag[i][j][k  ]=+2;
          iflag[i][j][k+1]=-2;
        }
      }
    }
  }

  /* ibody */
  for_aijk(i,j,k) {
    if(dom->ibody().off(i,j,k)) 
      iflag[i][j][k]=-1000;
  }

  /* trigger based on adens */
  for_ijk(i,j,k) {
    if(iflag[i][j][k] != -1000 && Interface(i,j,k)) {
      if(scp[i][j][k]<phisurf)
        iflag[i][j][k]=-1;
      else
        iflag[i][j][k]=+1; /* clr = 0.5? questionable...*/
    }
  }

  iflag.exchange();
}	

/***************************************************************************//**
*  Checks if the given cell is near an interface
*******************************************************************************/
bool VOF::Interface(const Sign dir, const Comp m,
                    const int i, const int j, const int k) {
  int of(0);
  if(dir==Sign::pos())
    of = 1;
 
  if(m==Comp::i())
    return boil::realistic(fs[m][i+of][j][k]);  
  else if(m==Comp::j())
    return boil::realistic(fs[m][i][j+of][k]);  
  else
    return boil::realistic(fs[m][i][j][k+of]);  

  return false;
}

bool VOF::Interface(const int i, const int j, const int k) {
  if(adens[i][j][k]>boil::pico) {
    return true;
  }

  return false;
}


