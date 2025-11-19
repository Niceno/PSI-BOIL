#include "momentum.h"

/******************************************************************************/
void Momentum::new_time_step(const Scalar * prs) {
  new_time_step(u,prs);
}

void Momentum::new_time_step(const Vector & vec, const Scalar * prs) {
  /*------------+
  |      dV  n  |
  |  f = -- u   |
  |      dt     |
  +------------*/
  for_m(m) 
    for_mijk(m,i,j,k) {
      const real rho = fluid()->rho(m,i,j,k);
      fold[m][i][j][k] = rho * dV(m,i,j,k) * vec[m][i][j][k] * time->dti();
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

#if 1
  /*------------+
  |        n-1  |
  |  f -= H     |
  |             |
  +------------*/
  if (conv_ts.Nm2()!=0.0)  {
    for_m(m) { 
      for_mijk(m,i,j,k) { /* conv_ts.Nm2() = -0.5 for adams-bashforth */
        fold[m][i][j][k] += conv_ts.Nm2() * cold[m][i][j][k]; 
      }
    }
  }

  /*------------+
  |       3  n  |
  |  f += - H   |
  |       2     |
  +------------*/
  /* a condition like: if(conv_ts != backward_euler()) would be good */
  if (conv_ts.Nm1()==0.0 && conv_ts.Nm2()==0.0) {
  } else {
    convection(&cold,prs);
    for_m(m) {
      for_mijk(m,i,j,k) /*conv_ts.Nm1()=1.5*/
        fold[m][i][j][k] += conv_ts.Nm1() * cold[m][i][j][k];
    }
  }
#endif
    
  /*------------+ 
  |       1  n  |
  |  f += - D   |
  |       2     |
  +------------*/
  /* a condition like: if(diff_ts != backward_euler()) would be good */
  if(diff_ts.Nm1() != 0.0)
    diffusion();
}
