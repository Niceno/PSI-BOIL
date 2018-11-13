#include "phasechange.h"
//#define ST_LENGTH
//#define DEBUG
using namespace std;

/******************************************************************************/
bool PhaseChange::incl_no_interface(const int i, const int j, const int k) {
/***************************************************************************//**
*  \brief check if cell(i,j,k) include interface
*  true: includes no interface,  false: includes interface
*******************************************************************************/
  real clrc = clr[i][j][k];
  real clrw = clr[i-1][j][k];
  real clre = clr[i+1][j][k];
  real clrs = clr[i][j-1][k];
  real clrn = clr[i][j+1][k];
  /* assume k-const is wall */
  //real clrb = clr[i][j][k-1];
  //real clrt = clr[i][j][k+1];

  /* west */
  if ((clrc-0.5)*(clrw-0.5)<0.0) {
    real legw = (0.5-clrc)/(clrw-clrc);
    if (legw <= 0.5) return false;
  }
  /* east */
  if ((clrc-0.5)*(clre-0.5)<0.0) {
    real lege = (0.5-clrc)/(clre-clrc);
    if (lege <= 0.5) return false;
  }
  /* south */
  if ((clrc-0.5)*(clrs-0.5)<0.0) {
    real legs = (0.5-clrc)/(clrs-clrc);
    if (legs <= 0.5) return false;
  }
  /* north */
  if ((clrc-0.5)*(clrn-0.5)<0.0) {
    real legn = (0.5-clrc)/(clrn-clrc);
    if (legn <= 0.5) return false;
  }

  return true;
}

