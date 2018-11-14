#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
real PhaseChangeVOF::gradtx(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in x-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  if(upwind_flag)
    return gradtx8(dir,i,j,k);
  else 
    return gradtx9(dir,i,j,k);
}

/******************************************************************************/
real PhaseChangeVOF::gradty(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in y-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  if(upwind_flag)
    return gradty8(dir,i,j,k);
  else 
    return gradty9(dir,i,j,k);
}

/******************************************************************************/
real PhaseChangeVOF::gradtz(const int dir, const int i, const int j, const int k){
/***************************************************************************//**
*  \brief calculate grad(tpr) in z-direction
*         dir= 1: grad(tpr) second order downwind
*         dir=-1: grad(tpr)  second order upwind
*******************************************************************************/
  if(upwind_flag)
    return gradtz8(dir,i,j,k);
  else 
    return gradtz9(dir,i,j,k);
}
