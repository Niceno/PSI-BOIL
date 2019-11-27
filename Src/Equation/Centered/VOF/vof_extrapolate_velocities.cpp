#include "vof.h"

//#define DEBUG

/******************************************************************************/
void VOF::extrapolate_velocities(const Scalar & scp, const Scalar & fext,
                                 const Matter * fluid, const Vector & umixed, 
                                 Vector & uliq, Vector & ugas,
                                 const ResRat & resrat) {
/***************************************************************************//**
*  \brief Calculate divergence-free liquid and gas velocities. Interface has to
*         be properly flagged.
*
*         input: clr, src, matter, uvw_mixed
*        output: uvw_liq, uvw_gas
*    throughput: fold -> mustn't be updated during the time step!
*     temporary: tempflag, stmp
*******************************************************************************/

  boil::timer.start("vof extrapolate velocities");
 
  /* step one: calculate tempflag */
  ev_flagging(scp,iflag,tempflag);

#ifdef DEBUG
  boil::plot->plot(scp, iflag, "c-iflag",time->current_step());
  boil::plot->plot(scp, tempflag, "c-pflag",time->current_step());
  //exit(0);
#endif

  /* step two: construct matrix. The unused "A" matrix is used */
  ev_discretize(fluid,tempflag,A);

#ifdef DEBUG
  boil::plot->plot(A.c,A.w,A.e,A.b,A.t, "ac-aw-ae-ab-at",time->current_step());
#endif

  /* step three: solve the system A*p = fext */
  const int niter = 200;
  //ev_solve(tempflag,A,fext,stmp,stmp2,false,niter);
  ev_solve(tempflag,A,fext,stmp,fold,true,niter,resrat);

  /* step four: project using the mixture-to-liquid velocity correction */
  for_m(m)
    for_avmijk(umixed,m,i,j,k)
      uliq[m][i][j][k] = umixed[m][i][j][k];
  ev_project(tempflag,fluid,stmp,uliq);

#ifdef DEBUG
  boil::plot->plot(uliq,scp, fext, stmp, "ucorr-c-f-p",time->current_step());
  for_ijk(i,j,k) {
    if(fabs(uliq.outflow(i,j,k)-umixed.outflow(i,j,k))>0.0)
      boil::oout<<i<<" "<<j<<" "<<k<<" | "
                <<uliq.outflow(i,j,k)*time->dti()<<" "<<umixed.outflow(i,j,k)*time->dti()
                <<" "<<fext[i][j][k]<<" | "
                <<uliq[Comp::i()][i][j][k]<<" "<<uliq[Comp::i()][i+1][j][k]<<" "
                <<uliq[Comp::k()][i][j][k]<<" "<<uliq[Comp::k()][i+1][j][k]<<" "
                <<boil::endl;
  }
  exit(0);
#endif

  boil::timer.stop("vof extrapolate velocities");

  return;
}

