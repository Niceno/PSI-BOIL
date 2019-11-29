#include "vof.h"

//#define DEBUG

/******************************************************************************/
void VOF::extrapolate_velocity(const Scalar & scp, const Scalar & fext,
                               const Matter * fluid, const Vector & umixed, 
                               Vector & unew, const ResRat & resrat,
                               const Sign & sig) {
/***************************************************************************//**
*  \brief Calculate divergence-free velocity field extension. Interface has to
*         be properly flagged.
*
*         sig = +1: obtain liquid velocity and vice versa
*
*         input: clr, src, matter, uvw_mixed
*        output: uvw_single_phase
*    throughput: pold_pos and pold_neg
*     temporary: tempflag, stmp
*******************************************************************************/

  boil::timer.start("vof extrapolate velocities");

#ifdef DEBUG
  interfacial_flagging(scp);
#endif

  /* step zero: choose direction */
  Scalar * pold;
  if(sig==Sign::neg()) {
    pold = &pold_neg;
  } else {
    pold = &pold_pos;
  }
 
  /* step one: calculate tempflag */
  ev_flagging(scp,iflag,tempflag,sig);

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
  const int niter = 2000;
  const bool converged = ev_solve(tempflag,A,fext,stmp,stmp2,false,niter,resrat);
  //const bool converged = ev_solve(tempflag,A,fext,stmp,*pold,true,niter,resrat);

  /* step four: project using the mixture-to-single velocity correction */
  for_m(m)
    for_avmijk(umixed,m,i,j,k)
      unew[m][i][j][k] = umixed[m][i][j][k];
  if(converged)
    ev_project(tempflag,fluid,stmp,unew);

#ifdef DEBUG
  boil::plot->plot(unew,scp, fext, stmp, "ucorr-c-f-p",time->current_step());
  for_ijk(i,j,k) {
    if(fabs(unew.outflow(i,j,k)-umixed.outflow(i,j,k))>0.0)
      boil::oout<<i<<" "<<j<<" "<<k<<" | "
                <<unew.outflow(i,j,k)*time->dti()<<" "<<umixed.outflow(i,j,k)*time->dti()
                <<" "<<fext[i][j][k]<<" | "
                <<unew[Comp::i()][i][j][k]<<" "<<unew[Comp::i()][i+1][j][k]<<" "
                <<unew[Comp::k()][i][j][k]<<" "<<unew[Comp::k()][i+1][j][k]<<" "
                <<boil::endl;
  }
  exit(0);
#endif

  boil::timer.stop("vof extrapolate velocities");

  return;
}

