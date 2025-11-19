#include "marching_cubes.h"

#define MG_LEVEL 3
/* 1 = nodes
   2 = nodes + faces
   3 = nodes + faces + edges
 */

/******************************************************************************/
void MarchingCubes::marker_gradient(Scalar & mga, const bool evalflag) {
/***************************************************************************//**
*  \brief Calculate the area density using marker gradient method.
*******************************************************************************/

#if MG_LEVEL == 1
  const real mult_node = 0.25; /* 4 per face */
  const real mult_face = 0.;
  const real mult_edge = 0.;
#elif MG_LEVEL == 2
  const real mult_node = 0.5*0.25; /* face takes 1/2 */
  const real mult_face = 0.5;
  const real mult_edge = 0.;
#elif MG_LEVEL == 3
  const real mult_node = 0.25*0.25; /* face takes 1/4 and edges 1/2 */
  const real mult_face = 0.25;
  const real mult_edge = 0.5*0.25; /* 4 per face */
#endif

  if(evalflag)
    evaluate_nodes();

  for_vijk(mga,i,j,k) {
    real W = (*clr)[i][j][k]+(*clr)[i-1][j][k] >= 2.*clrsurf;
    real E = (*clr)[i][j][k]+(*clr)[i+1][j][k] >= 2.*clrsurf;
    real S = (*clr)[i][j][k]+(*clr)[i][j-1][k] >= 2.*clrsurf;
    real N = (*clr)[i][j][k]+(*clr)[i][j+1][k] >= 2.*clrsurf;
    real B = (*clr)[i][j][k]+(*clr)[i][j][k-1] >= 2.*clrsurf;
    real T = (*clr)[i][j][k]+(*clr)[i][j][k+1] >= 2.*clrsurf;

    real WS = (*clr)[i][j][k]+(*clr)[i-1][j-1][k]+(*clr)[i-1][j][k]+(*clr)[i][j-1][k] >= 4.*clrsurf;
    real WN = (*clr)[i][j][k]+(*clr)[i-1][j+1][k]+(*clr)[i-1][j][k]+(*clr)[i][j+1][k] >= 4.*clrsurf;
    real ES = (*clr)[i][j][k]+(*clr)[i+1][j-1][k]+(*clr)[i+1][j][k]+(*clr)[i][j-1][k] >= 4.*clrsurf;
    real EN = (*clr)[i][j][k]+(*clr)[i+1][j+1][k]+(*clr)[i+1][j][k]+(*clr)[i][j+1][k] >= 4.*clrsurf;
    real WB = (*clr)[i][j][k]+(*clr)[i-1][j][k-1]+(*clr)[i-1][j][k]+(*clr)[i][j][k-1] >= 4.*clrsurf;
    real WT = (*clr)[i][j][k]+(*clr)[i-1][j][k+1]+(*clr)[i-1][j][k]+(*clr)[i][j][k+1] >= 4.*clrsurf;
    real EB = (*clr)[i][j][k]+(*clr)[i+1][j][k-1]+(*clr)[i+1][j][k]+(*clr)[i][j][k-1] >= 4.*clrsurf;
    real ET = (*clr)[i][j][k]+(*clr)[i+1][j][k+1]+(*clr)[i+1][j][k]+(*clr)[i][j][k+1] >= 4.*clrsurf;
    real SB = (*clr)[i][j][k]+(*clr)[i][j-1][k-1]+(*clr)[i][j-1][k]+(*clr)[i][j][k-1] >= 4.*clrsurf;
    real ST = (*clr)[i][j][k]+(*clr)[i][j-1][k+1]+(*clr)[i][j-1][k]+(*clr)[i][j][k+1] >= 4.*clrsurf;
    real NB = (*clr)[i][j][k]+(*clr)[i][j+1][k-1]+(*clr)[i][j+1][k]+(*clr)[i][j][k-1] >= 4.*clrsurf;
    real NT = (*clr)[i][j][k]+(*clr)[i][j+1][k+1]+(*clr)[i][j+1][k]+(*clr)[i][j][k+1] >= 4.*clrsurf;

    real EST = nodalvals[i+1][j  ][k+1] >= clrsurf;
    real WST = nodalvals[i  ][j  ][k+1] >= clrsurf;
    real ESB = nodalvals[i+1][j  ][k  ] >= clrsurf;
    real WSB = nodalvals[i  ][j  ][k  ] >= clrsurf;
    real ENT = nodalvals[i+1][j+1][k+1] >= clrsurf;
    real WNT = nodalvals[i  ][j+1][k+1] >= clrsurf;
    real ENB = nodalvals[i+1][j+1][k  ] >= clrsurf;
    real WNB = nodalvals[i  ][j+1][k  ] >= clrsurf;

    mga[i][j][k] = sqrt( pow((mult_face*(E-W)
                        +mult_node*(ENT+ENB+EST+ESB
                                -WNT-WNB-WST-WSB)
                        +mult_edge*(ES+EN+EB+ET-WS-WN-WB-WT)
                        )/mga.dxc(i),2.0)
                       + pow((mult_face*(N-S)
                        +mult_node*(ENT+ENB+WNT+WNB
                                -EST-ESB-WST-WSB)
                        +mult_edge*(WN+EN+NB+NT-WS-ES-SB-ST)
                        )/mga.dyc(j),2.0)
                       + pow((mult_face*(T-B)
                        +mult_node*(EST+WST+ENT+WNT
                                -ESB-WSB-ENB-WNB)
                        +mult_edge*(WT+ET+ST+NT-WB-EB-SB-NB)
                        )/mga.dzc(k),2.0) );
  }
  mga.bnd_update();
  mga.exchange_all();

  return;
}
