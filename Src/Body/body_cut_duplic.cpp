#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"

/******************************************************************************/
int Body::cut_duplic(Plane * p,
                     const real x[], const real y[], const real z[],
                     const int n_in_fluid[][2][2]) const {
/***************************************************************************//**
*  \brief check duplication of cell-face and cut-plane.
*         return value: 0 = duplication, cell is solid
*                       1 = no duplication
*                       2 = duplication, cell is fluid
*******************************************************************************/

  /* cell center */
  real xc = 0.5*(x[0]+x[1]);
  real yc = 0.5*(y[0]+y[1]);
  real zc = 0.5*(z[0]+z[1]);

  /* i-direction */
  for(int i=0; i<2; i++){
    int isum=0;
    for(int j=0; j<2; j++){
      for(int k=0; k<2; k++){
        isum += n_in_fluid[i][j][k];
      }
    }
    if(isum==8){  // duplication
      if( p->distance(xc,yc,zc) < 0.0){ // solid cell
        return 0;
      } else {
        return 2;
      }
    }
  }

  /* y-direction */
  for(int j=0; j<2; j++){
    int isum=0;
    for(int i=0; i<2; i++){
      for(int k=0; k<2; k++){
        isum += n_in_fluid[i][j][k];
      }
    }
    if(isum==8){  // duplication
      //std::cout<<"cut_duplic: j-duplication\n";
      if( p->distance(xc,yc,zc) < 0.0){ // solid cell
        return 0;
      } else {
        return 2;
      }
    }
  }

  /* z-direction */
  for(int k=0; k<2; k++){
    int isum=0;
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        isum += n_in_fluid[i][j][k];
      }
    }
    if(isum==8){  // duplication
      //std::cout<<"cut_duplic: k-duplication\n";
      if( p->distance(xc,yc,zc) < 0.0){ // solid cell
        return 0;
      } else {
        return 2;
      }
    }
  }
  return 1;
}

/*-----------------------------------------------------------------------------+
 '$Id: body_cut_duplic.cpp,v 1.1 2011/03/28 08:12:04 sato Exp $'/
+-----------------------------------------------------------------------------*/
