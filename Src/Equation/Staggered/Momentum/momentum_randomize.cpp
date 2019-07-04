#include "momentum.h"                

/******************************************************************************/
void Momentum::randomize(RandomFlow & rf, const real magn) {

  /*-------------------------------+
  |  create random velocity field  |
  +-------------------------------*/
  real max=0.0;
  real x[3];
  for_m(m) {
    if(   m==Comp::u() && ifull
       || m==Comp::v() && jfull
       || m==Comp::w() && kfull
      )
      for_vmijk(u,m,i,j,k) {
        x[0] = u.xc(m,i); 
        x[1] = u.yc(m,j); 
        x[2] = u.zc(m,k);
        rf.get_vector( time.current_time(), x, &u[m][i][j][k], m );
        max = maxr(u[m][i][j][k], max);
      }
  }

  /*-------------------------+ 
  |  scale the fluctuations  |
  +-------------------------*/
  real factor = magn/max;
  for_m(m) {
    for_vmijk(u,m,i,j,k) {
      u[m][i][j][k] *= factor;
    }
  }
}
