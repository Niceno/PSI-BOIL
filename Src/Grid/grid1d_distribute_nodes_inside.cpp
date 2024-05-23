#include "grid1d.h"

/******************************************************************************/
void Grid1D::distribute_nodes_inside(const real & x1, const real & xn,
                                     const real & D1, const real & DN) {

  int    nc_in_m = nc_in-1;
  real L  = xn - x1;

  assert(nc_in>1);

/*-----------------------+        nc_in = 2
|                        |
|  special case: nc_in == 2  |    0   1   2   3
|                        |  | o |-O-+-O-| o |
+------------------------+  0   1   2   3   4  */
  if(nc_in == 2) {
    real dx = L*0.5;

    x_node[boil::BW]   = x1;
    x_node[boil::BW+2] = xn;
    x_node[boil::BW+1] = 0.5 * (x1 + xn);
    return;
  }

/*---------------------------------------------------------------+
|                                                                |
|  Grid distribution is determined from the following function   |
|                                                                |
|  x = a + b y + c y^2 + d y^3                                   |
|                                                                |
|  One should imagine this function as taking integer arguments. |
|                                                                |
|  conditions:                                                   |
|  ~~~~~~~~~~~                                                   |
|                                                                |
|  1. x(0) = 0                                                   |
|                                                                |
|     a = 0                                                      |
|                                                                |
|  2. x(1) = D1                                                  |
|                                                                |
|     b + c + d = D1                                             |
|                                                                |
|  3. X(nc_in) = L                                               |
|                                                                |
|     b nc_in + c nc_in^2 + d nc_in^3 = L                        |
|                                                                |
|  4. X(nc_in-1) = L - DN 
|                                                                |
|     b (nc_in-1) + c (nc_in-1)^2 + d (nc_in-1)^3 = L - DN       |
|                                                                |
|  2,3,4 =>                                                      |
|                                                                |
|  |1        1           1         | |b|   |D1  |                |
|  |nc_in    nc_in^2     nc_in^3   | |c| = |L   |                |
|  |nc_in-1 (nc_in-1)^2 (nc_in-1)^3| |d|   |L-DN|                |
|                                                                |
+---------------------------------------------------------------*/
  const int n = 3;
  real A[n][n];
  real bcd[n];
  real f[n];

  A[0][0]=1;        A[0][1]=1.0;              A[0][2]=1.0;
  A[1][0]=nc_in;    A[1][1]=nc_in*nc_in;      A[1][2]=nc_in*nc_in*nc_in;
  A[2][0]=nc_in_m;  A[2][1]=nc_in_m*nc_in_m;  A[2][2]=nc_in_m*nc_in_m*nc_in_m;

  f[0] = D1;
  f[1] = L; 
  f[2] = L-DN;

  /*--------------------+ 
  |                     |
  |  Gauss elimination  |
  |                     |
  +--------------------*/

  /*----------+
  |  forward  |
  +----------*/
  for(int i=0; i<n; i++) {

    /* elimination of bellow diagonal */
    for(int j=i+1; j<n; j++) {
      real r = A[j][i]/A[i][i];

      /* handle row */
      for(int k=0; k<n; k++) 
        A[j][k] -= r * A[i][k];

      f[j] -= r * f[i];
    }
  }

  /*-----------+
  |  backward  |
  +-----------*/
  bcd[n-1] = f[n-1]/A[n-1][n-1];

  for(int i=n-2; i>=0; i--) {
    real v = f[i];

    for(int j=i+1; j<n; j++)
      v -= A[i][j]*bcd[j];

    bcd[i] = v/A[i][i];
  }

  /*--------------------------+            nc_in = 4
  |                           |
  |  create node coordinates  |    0   1   2   3   4   5
  |                           |  | o |-O-+-O-+-O-+-O-| o |
  +---------------------------+  0   1   2   3   4   5   6  */
  for(int i=0; i<nc_in+1; i++) 
    x_node[i+boil::BW] = x1 + bcd[0] * i + bcd[1] * (i*i) + bcd[2] * (i*i*i);

  /*-----------------------------+
  |  check if they are monotone  |
  +-----------------------------*/
  for(int i=boil::BW+1; i<nc_in+2; i++) 
    if(x_node[i] <= x_node[i-1]) {
      boil::aout << "Failure! Probably due to too big D1 and DN" << boil::endl;
      exit(0);
    }
}
