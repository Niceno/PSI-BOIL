#include "momentum.h"

/******************************************************************************/
void Momentum::new_time_step() {

  /*------------+
  |      dV  n  |
  |  f = -- u   |
  |      dt     |
  +------------*/
  for_m(m) 
    for_mijk(m,i,j,k) {
      const real rho = fluid()->rho(m,i,j,k);
      fold[m][i][j][k] = rho * dV(m,i,j,k) * u[m][i][j][k] * time->dti();
    }

  /* correct for volumes in immersed boundary */
  for_m(m) {
    if(dom->ibody().nccells(m) > 0) {
      for(int cc=0; cc<dom->ibody().nccells(m); cc++)
        if( dom->ibody().cut(m,cc) ) {
          int i,j,k;
          dom->ibody().ijk(m,cc,&i,&j,&k); 
          fold[m][i][j][k] *= dom->ibody().fV(m,cc); 
        }
    }
  }

  /*------------+
  |        n-1  |
  |  f -= H     |
  |             |
  +------------*/
  for_m(m) 
    for_mijk(m,i,j,k) {
                       /* conv_ts.Nm2() = -0.5 for adams-bashforth */
      fold[m][i][j][k] += conv_ts.Nm2() * cold[m][i][j][k]; 
    }

  /*------------+
  |       3  n  |
  |  f += - H   |
  |       2     |
  +------------*/
  /* a condition like: if(conv_ts != backward_euler()) would be good */
  convection(&cold);
  for_m(m)
    for_mijk(m,i,j,k)
      fold[m][i][j][k] += conv_ts.Nm1() * cold[m][i][j][k]; /*conv_ts.Nm1()=1.5*/
    
  /*------------+ 
  |       1  n  |
  |  f += - D   |
  |       2     |
  +------------*/
  /* a condition like: if(diff_ts != backward_euler()) would be good */
  diffusion();
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_new_time_step.cpp,v 1.17 2010/03/25 08:15:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
