#include "momentum.h"
#include "../../../Plot/plot.h"

/******************************************************************************/
void Momentum::solve(const ResRat & factor) {

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
        }
      }
    }
  } /* is there an immersed body */

  /*--------+
  |  solve  |
  +--------*/
  for_m(m) {

    Matrix * Am = A[~m];

    if(m==Comp::u()) 
      solver->solve(*Am, u(m), fnew(m), 
                     MaxIter(10), "u", factor);
    if(m==Comp::v()) 
      solver->solve(*Am, u(m), fnew(m), 
                     MaxIter(10), "v", factor);
    if(m==Comp::w()) 
      solver->solve(*Am, u(m), fnew(m), 
                     MaxIter(10), "w", factor);
  }

  /* set velocity in solid zero */
  if(dom->ibody().nccells() > 0) {
    for_m(m){
      for_mijk(m,i,j,k)
        if(dom->ibody().off(m,i,j,k))u[m][i][j][k]=0.0;
    }
  }

  insert_bc();
  u.exchange_all();

  boil::timer.stop("momentum solver");

  // boil::plot->plot(fnew, "uvw-fnew", time->current_step());
  // boil::plot->plot(cnew, "uvw-cnew", time->current_step());

  /*---------------------------------+
  |  scale velocities at the outlet  |
  +---------------------------------*/
  outlet();
}
