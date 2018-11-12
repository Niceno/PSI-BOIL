#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::globals(int * i, int * j, int * k) const {
/*---------------------------------------+
|  sends local i,j,k but returns global  |
+---------------------------------------*/

  *i = *i+cr_x.first()-1;
  *j = *j+cr_y.first()-1;
  *k = *k+cr_z.first()-1;
}

/******************************************************************************/
int Domain::global_I(int i) const {
/*---------------------------------------+
|  receives local i, but returns global  |
+---------------------------------------*/

  return i+cr_x.first()-1;
}

/******************************************************************************/
int Domain::global_J(int j) const {
/*---------------------------------------+
|  receives local j, but returns global  |
+---------------------------------------*/

  return j+cr_y.first()-1;
}

/******************************************************************************/
int Domain::global_K(int k) const {
/*---------------------------------------+
|  receives local k, but returns global  |
+---------------------------------------*/

  return k+cr_z.first()-1;
}
