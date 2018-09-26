#include "scalar.h"

/******************************************************************************/
real Scalar::average_I(int I) const {
/*------------------------------------------+
|  capital "I" means it is a global number  |
+------------------------------------------*/
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_I(I) ) {

    int i=dom->local_i(I);

    for_jk(j,k) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_I(I)) ) return 0.0;
 
  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_J(int J) const {
/*------------------------------------------+
|  capital "J" means it is a global number  |
+------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_J(J) ) {

    int j=dom->local_j(J);

    for_ik(i,k) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_J(J)) ) return 0.0;
 
  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_K(int K) const {
/*------------------------------------------+
|  capital "K" means it is a global number  |
+------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_K(K) ) {

    int k=dom->local_k(K); // only i will be needed

    for_ij(i,j) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_K(K)) ) return 0.0;
 
  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_IJ(int I, int J) const {
/*-----------------------------------------------------+
|  capital "I" and "J" mean it is a global coordinate  |
+-----------------------------------------------------*/
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_I(I) && dom->contains_J(J) ) {

    int i=dom->local_i(I);
    int j=dom->local_j(J);

    for_k(k) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_I(I) && dom->contains_J(J)) ) return 0.0;

  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_IK(int I, int K) const {
/*-----------------------------------------------------+
|  capital "I" and "K" mean it is a global coordinate  |
+-----------------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_I(I) && dom->contains_K(K) ) {

    int i=dom->local_i(I);
    int k=dom->local_k(K);

    for_j(j) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_I(I) && dom->contains_K(K)) ) return 0.0;

  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_JK(int J, int K) const {
/*-----------------------------------------------------+
|  capital "J" and "K" mean it is a global coordinate  |
+-----------------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  if( dom->contains_J(J) && dom->contains_K(K) ) {

    int j=dom->local_j(J); // only i will be needed
    int k=dom->local_k(K); // only i will be needed

    for_i(i) {
      sum += val[i][j][k];
      count++;
    } 
  }

  boil::cart.sum_real(&sum);
  boil::cart.sum_real(&count);

  if( !(dom->contains_J(J) && dom->contains_K(K)) ) return 0.0;

  assert(count != 0.0);

  return sum/count;
}

/******************************************************************************/
real Scalar::average_i(int i) const {
/*------------------------------------------+
|  low-case "i" means it is a local number  |
+------------------------------------------*/
  real count = 0.0;
  real sum   = 0.0;

  for_jk(j,k) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/******************************************************************************/
real Scalar::average_j(int j) const {
/*------------------------------------------+
|  low-case "j" means it is a local number  |
+------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  for_ik(i,k) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/******************************************************************************/
real Scalar::average_k(int k) const {
/*------------------------------------------+
|  low-case "k" means it is a local number  |
+------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  for_ij(i,j) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/******************************************************************************/
real Scalar::average_ij(int i, int j) const {
/*-----------------------------------------------------+
|  low-case "i" and "j" mean it is a local coordinate  |
+-----------------------------------------------------*/
  real count = 0.0;
  real sum   = 0.0;

  for_k(k) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/******************************************************************************/
real Scalar::average_ik(int i, int k) const {
/*-----------------------------------------------------+
|  low-case "i" and "k" mean it is a local coordinate  |
+-----------------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  for_j(j) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/******************************************************************************/
real Scalar::average_jk(int j, int k) const {
/*-----------------------------------------------------+
|  low-case "i" and "k" mean it is a local coordinate  |
+-----------------------------------------------------*/
  
  real count = 0.0;
  real sum   = 0.0;

  for_i(i) {
    sum += val[i][j][k];
    count++;
  } 

  return sum/count;
}

/*-----------------------------------------------------------------------------+
 '$Id: scalar_average.cpp,v 1.5 2012/03/13 09:24:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
