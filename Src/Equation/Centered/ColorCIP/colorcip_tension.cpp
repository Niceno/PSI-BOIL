#include "colorcip.h"
#include <iomanip>

/******************************************************************************/
void ColorCIP::tension(Vector * vec, const Matter matt) {
/***************************************************************************//**
*  \brief Calculate surface tension
*         Algorithm
*           1st step: calculate smoothed color function
*           2nd step: calculate curvature
*           3rd step: calculate body force
*         Variables
*           color function          : clr
*           smoothed color function : clrn (use only for curvature calculation)
*           curvature               : kappa
*           body force              : vec
*******************************************************************************/

  /*-----------+
  |  1st step  |
  +-----------*/
#if 1
  smooth(clr, clrn, 5);
#else
  smooth_ls(clr, clrn, 16);
#endif
  insert_bc(clrn);
  clrn.exchange_all();

  /*-----------+
  |  2nd step  |
  +-----------*/
  /* calculate curvature */
  curv(clrn);

  /* wall adhesion */
  const real cangle = 45;
  theta = cangle/ 180.0 * pi;
  bdcurv(clrn,theta);
  kappa.exchange();

#if 0
  for_ijk(i,j,k)    // convert distance function to color function
                    // for the case smooth_ls is used above.
    clrn[i][j][k]=0.5+0.5*tanh(clrn[i][j][k]/(sqrt(2.0)*ww));
#endif

  /*-----------+
  |  3rd step  |
  +-----------*/
  /* calculate body force */
  Comp m;
  m = Comp::u();
  for_vmijk((*vec),m,i,j,k) {
    (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                       * 0.5*(kappa[i-1][j][k]+kappa[i][j][k])
                       * (clrn[i][j][k]-clrn[i-1][j][k])/vec->dxc(m,i)
                       * vec->dV(m,i,j,k);
  }

  m = Comp::v();
  for_vmijk((*vec),m,i,j,k) {
    (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                       * 0.5*(kappa[i][j-1][k]+kappa[i][j][k])
                       * (clrn[i][j][k] - clrn[i][j-1][k])/vec->dxc(m,i)
                       * vec->dV(m,i,j,k);
  }

  m = Comp::w();
  for_vmijk((*vec),m,i,j,k) {
    (*vec)[m][i][j][k] += matt.sigma(m,i,j,k)
                       * 0.5*(kappa[i][j][k-1]+kappa[i][j][k])
                       * (clrn[i][j][k] - clrn[i][j][k-1])/vec->dxc(m,i)
                       * vec->dV(m,i,j,k);
  }
  vec->exchange();
}
