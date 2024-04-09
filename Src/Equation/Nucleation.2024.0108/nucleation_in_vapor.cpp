#include "nucleation.h"

/******************************************************************************/
bool Nucleation::in_vapor(const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  //if(matter_sig==Sign::pos()) {
    //return cht->topo->below_interface(i,j,k);
    return (*clr)[i][j][k]<threshold_c;
  //} else {
  //  return cht->topo->above_interface(i,j,k);
  //}
}

/******************************************************************************/
bool Nucleation::in_vapor(const real c) const {
/***************************************************************************//**
*  \brief test if in vapor taking into account sign
*******************************************************************************/
  //if(matter_sig==Sign::pos()) {
    //return cht->topo->below_interface(c);
    return c<threshold_c;
  //} else {
  //  return cht->topo->above_interface(c);
  //}
}

/******************************************************************************/
bool Nucleation::below_threshold(const real c) const {
/***************************************************************************//**
*  \brief test if c below the threshold
*******************************************************************************/
  //if(matter_sig==Sign::pos()) {
    return c<threshold_c;
  //} else {
  //  return c>threshold_c;
  //}
}

/******************************************************************************/
bool Nucleation::below_threshold(const int i, const int j, const int k) const {
/***************************************************************************//**
*  \brief test if color below threshold
*******************************************************************************/
  return below_threshold((*clr)[i][j][k]);
}
