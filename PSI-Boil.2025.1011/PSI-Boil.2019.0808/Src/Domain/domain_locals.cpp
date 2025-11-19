#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::locals(int * i, int * j, int * k) const {
/*---------------------------------------+
|  sends global i,j,k but returns local  |
+---------------------------------------*/

  if( contains_IJK(*i, *j, *k) ) {
    *i = *i-cr_x.first()+1;
    *j = *j-cr_y.first()+1;
    *k = *k-cr_z.first()+1;
  }
  else {
    *i = -1;
    *j = -1;
    *k = -1;
  }
}

/******************************************************************************/
int Domain::local_i(int I) const {
/*---------------------------------------+
|  receives global i, but returns local  |
+---------------------------------------*/

  if( contains_I(I) ) return I-cr_x.first()+1;
  else                return -1;
}

/******************************************************************************/
int Domain::local_j(int J) const {
/*---------------------------------------+
|  receives global j, but returns local  |
+---------------------------------------*/

  if( contains_J(J) ) return J-cr_y.first()+1;
  else                return -1;
}

/******************************************************************************/
int Domain::local_k(int K) const {
/*---------------------------------------+
|  receives global k, but returns local  |
+---------------------------------------*/

  if( contains_K(K) ) return K-cr_z.first()+1;
  else                return -1;
}
