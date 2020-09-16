#include "vof.h"

real VOF::extract_cl_velocity_2d(const Comp ctangential, const Comp cnormal,
                                 const Sign sig, 
                                 int * IG, int * PN,
                                 const Range<int> ridx) {
/***************************************************************************//**
*  \brief Extract contact line velocity for a 2D system. 
*
*         Limitations: only single bubble/droplet with one contact line
*                      is assumed (= no significant inward bending of 
*                      interface), wall assumed in negative z-dir.
*
*******************************************************************************/

  assert(ctangential==Comp::i());
  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());

  bool ranged(false);
  if(ridx.exists()) {
    ranged = true;
  }

  real veloc(boil::unreal);
  int Iglob(boil::unint);
  int procnum(boil::unint);

  for_i(i) {
    int Iloc = phi.domain()->global_I(i)-boil::BW+1;
    if(Iloc<Iglob && (!ranged||ridx.contains(Iloc))) {
      for_jk(j,k) {
        if(dom->ibody().on(i,j,k)) {
          if(dom->ibody().off(i,j,k-1) || (k==sk() && kminw)) {
            if(abs(iflag[i][j][k])==1) {
              Iglob = Iloc;
              veloc = 0.5*((*u)[Comp::u()][i][j][k]+(*u)[Comp::u()][i+1][j][k]);
            } /* at the interface */
          } /* next to wall */
        } /* ibody on */
      } /* jk */
    } /* Iloc smaller that Inglob, we are in range */
  } /* local i */

  /* figure out the processor with lowest Iglob */
  int Iglob_temp = Iglob;
  boil::cart.min_int(&Iglob);

  /* collect velocity and processor number */
  if(Iglob==Iglob_temp) {
    procnum = boil::cart.iam();
  } else {
    veloc = boil::unreal;
  }
  boil::cart.min_int(&procnum);
  boil::cart.min_real(&veloc);

  /* return */
  if(boil::realistic(Iglob)) {
    if(IG)
      *IG = Iglob;
    if(PN)
      *PN = procnum;
    return veloc;
  } else {
    if(IG)
      *IG = 0;
    if(PN)
      *PN = 0;
    return 0.;
  }
}
