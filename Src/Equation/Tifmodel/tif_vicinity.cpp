#include "tif.h"

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
bool TIF::Vicinity(const int i, const int j, const int k) {

  if(Interface(i+1,j,k)) {
    return true;
  } 
  if(Interface(i-1,j,k)) {
    return true;
  } 
  if(Interface(i,j+1,k)) {
    return true;
  } 
  if(Interface(i,j-1,k)) {
    return true;
  } 
  if(Interface(i,j,k+1)) {
    return true;
  } 
  if(Interface(i,j,k-1)) {
    return true;
  } 
  
  return false;
}
