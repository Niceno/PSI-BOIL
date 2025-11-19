#include "momentum.h"
#include "../../../Plot/plot.h"

/******************************************************************************/
void Momentum::solve(const ResTol & toler, const ResRat & factor) {
  solve_wo_outlet(toler,factor);

  /*---------------------------------+
  |  scale velocities at the outlet  |
  +---------------------------------*/
  outlet();

  return;
}

/******************************************************************************/
void Momentum::solve_wo_outlet(const ResTol & toler, const ResRat & factor) {

  boil::timer.start("momentum solver");

  /*-----------------------------+
  |  create final source vector  |
  +-----------------------------*/
  for_m(m) {
    int ip(0), jp(0), kp(0);
    if(m==Comp::u()) ip=-1;
    if(m==Comp::v()) jp=-1;
    if(m==Comp::w()) kp=-1;
    
    for_mijk(m,i,j,k) {
      fnew[m][i][j][k] = fold[m][i][j][k]
                       + cnew[m][i][j][k] * conv_ts.N()
                       + gradp[m][i][j][k] 
                       + fext[m][i][j][k];
    }
  }

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) {
    if(!ib_trust_vel_wall) {
    /* a very crude cut-off */
    for_m(m) {
      int ip=0, jp=0, kp=0;
      if(m==Comp::u()) ip++;
      if(m==Comp::v()) jp++;
      if(m==Comp::w()) kp++;
      for_mijk(m,i,j,k) {
        if( dom->ibody().off(i,   j,   k   ) || 
            dom->ibody().off(i-ip,j-jp,k-kp) ) {
          fnew[m][i][j][k] = 0.0;
          //fnew[m][i][j][k] = u[m][i][j][k];
        }
      }
    }
    }
  } /* is there an immersed body */

  /*--------+
  |  solve  |
  +--------*/
  for_m(m) {

    Matrix * Am = A[~m];

    if(m==Comp::u()) {
      if(ifull) {
        solver->solve(*Am, u(m), fnew(m), min_iter, 
                       MaxIter(10), "u", factor,toler,scale*time->dti());
      } else {
        for_avmijk(u,m,i,j,k) 
          u[m][i][j][k] = 0.0;
      }
    }
    if(m==Comp::v()) { 
      if(jfull) {
        solver->solve(*Am, u(m), fnew(m), min_iter,
                       MaxIter(10), "v", factor,toler,scale*time->dti());
      } else {
        for_avmijk(u,m,i,j,k) 
          u[m][i][j][k] = 0.0;
      }
    }
    if(m==Comp::w()) {
      if(kfull) {
        solver->solve(*Am, u(m), fnew(m), min_iter,
                       MaxIter(10), "w", factor,toler,scale*time->dti());
      } else {
        for_avmijk(u,m,i,j,k) 
          u[m][i][j][k] = 0.0;
      }
    }
  }

  /* set velocity in solid zero */
  if(!ib_trust_vel_wall) {
    if(dom->ibody().nccells() > 0) {
      for_m(m){
        for_mijk(m,i,j,k)
          if(dom->ibody().off(m,i,j,k)) {
            u[m][i][j][k]=0.0;
          }
      }
    }
  }

  u.bnd_update_nooutlet();
  u.exchange_all();

  boil::timer.stop("momentum solver");

  // boil::plot->plot(fnew, "uvw-fnew", time->current_step());
  // boil::plot->plot(cnew, "uvw-cnew", time->current_step());

  return;
}
