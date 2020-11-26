#include "global_vectorproduct.h"

namespace boil {

/******************************************************************************/
void crossProduct(real vect_C[], real vect_A[], real vect_B[]){
  vect_C[0] = vect_A[1]*vect_B[2] - vect_A[2]*vect_B[1];
  vect_C[1] = vect_A[2]*vect_B[0] - vect_A[0]*vect_B[2];
  vect_C[2] = vect_A[0]*vect_B[1] - vect_A[1]*vect_B[0];
}

/******************************************************************************/
real dotProduct(real vect_A[], real vect_B[]){
  real product = 0; 
  for (int i = 0; i < 3; i++){ 
    product = product + vect_A[i] * vect_B[i];
  }
  return product; 
} 

}

