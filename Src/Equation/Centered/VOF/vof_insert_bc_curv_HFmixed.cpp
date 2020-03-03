#include "vof.h"

//#define DEBUG

/******************************************************************************/
void VOF::insert_bc_curv_HFmixed(const Scalar & scp,
                                 const Comp ctangential, const Comp cnormal,
                                 const Sign sig) {
/***************************************************************************//**
*  \brief Calculate curvature using hybrid height-function/divergence-of-normal
*         approach in 2D geometry. The detachment criterion is set.
*
*         Reference: derived by me (Lubomir)
*
*         Limitations: wall assumed in negative z-dir.
*
*     temporary: tempflag2
*     output: kappa
*******************************************************************************/
  
  tempflag2 = -1000;

  assert(ctangential==Comp::i());
  assert(cnormal==Comp::k());
  assert(sig==Sign::neg());


  /**** step 1: flag bottom-layer cells ****/
  real scpref_aligned = mult_wall < 0 ? tol_wall : 1.-tol_wall;
  real scpref_nonaligned = mult_wall < 0 ? 1.-tol_wall : tol_wall;

#ifdef DEBUG
  boil::oout<<mult_wall<<" "<<scpref_aligned<<" "<<scpref_nonaligned<<boil::endl;
#endif

  int aligned_0_min, aligned_1_min, aligned_2_min;
  int nonaligned_0_min, nonaligned_1_min, nonaligned_2_min;
  int nontrivial_0_min, nontrivial_1_min, nontrivial_2_min;

  /* boil::unint is the integer version of boil::unreal */
  aligned_0_min = aligned_1_min = aligned_2_min          = boil::unint;
  nonaligned_0_min = nonaligned_1_min = nonaligned_2_min = boil::unint;
  nontrivial_0_min = nontrivial_1_min = nontrivial_2_min = boil::unint;

  int aligned_0_max, aligned_1_max, aligned_2_max;
  int nonaligned_0_max, nonaligned_1_max, nonaligned_2_max;
  int nontrivial_0_max, nontrivial_1_max, nontrivial_2_max;

  aligned_0_max = aligned_1_max = aligned_2_max          = -boil::unint;
  nonaligned_0_max = nonaligned_1_max = nonaligned_2_max = -boil::unint;
  nontrivial_0_max = nontrivial_1_max = nontrivial_2_max = -boil::unint;
  
  /*
     1: color of cell corresponds to the one expected within tolerance
    -1: color of cell corresponds to the inverse phase
     0: non-trivial value of color
  */

  int * aligned_min, * aligned_max;
  int * nonaligned_min, * nonaligned_max;
  int * nontrivial_min, * nontrivial_max;

  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      bool flag(false);
      if       (dom->ibody().off(i,j,k-1)||(kminw && (k==sk()  ))) {
        flag = true;
        aligned_min = &aligned_0_min;
        aligned_max = &aligned_0_max;
        nonaligned_min = &nonaligned_0_min;
        nonaligned_max = &nonaligned_0_max;
        nontrivial_min = &nontrivial_0_min;
        nontrivial_max = &nontrivial_0_max;
      } else if(dom->ibody().off(i,j,k-2)||(kminw && (k==sk()+1))) {
        flag = true;
        aligned_min = &aligned_1_min;
        aligned_max = &aligned_1_max;
        nonaligned_min = &nonaligned_1_min;
        nonaligned_max = &nonaligned_1_max;
        nontrivial_min = &nontrivial_1_min;
        nontrivial_max = &nontrivial_1_max;
      } else if(dom->ibody().off(i,j,k-3)||(kminw && (k==sk()+2))) {
        flag = true;
        aligned_min = &aligned_2_min;
        aligned_max = &aligned_2_max;
        nonaligned_min = &nonaligned_2_min;
        nonaligned_max = &nonaligned_2_max;
        nontrivial_min = &nontrivial_2_min;
        nontrivial_max = &nontrivial_2_max;
      }
      if(flag) {
        if       ( mult_wall*(scp[i][j][k]-scpref_aligned   )>=0.) {
          tempflag2[i][j][k] =  1;
          *aligned_min = std::min(*aligned_min,scp.domain()->global_I(i));
          *aligned_max = std::max(*aligned_max,scp.domain()->global_I(i));
        } else if(-mult_wall*(scp[i][j][k]-scpref_nonaligned)>=0.) {
          tempflag2[i][j][k] = -1;
          *nonaligned_min = std::min(*nonaligned_min,scp.domain()->global_I(i));
          *nonaligned_max = std::max(*nonaligned_max,scp.domain()->global_I(i));
        } else {
          tempflag2[i][j][k] =  0;
          *nontrivial_min = std::min(*nontrivial_min,scp.domain()->global_I(i));
          *nontrivial_max = std::max(*nontrivial_max,scp.domain()->global_I(i));
        }
      }
    }
  }
  tempflag2.bnd_update_symmetry();
  tempflag2.exchange();

  boil::cart.min_int(&aligned_0_min);
  boil::cart.min_int(&aligned_1_min);
  boil::cart.min_int(&aligned_2_min);

  boil::cart.min_int(&nonaligned_0_min);
  boil::cart.min_int(&nonaligned_1_min);
  boil::cart.min_int(&nonaligned_2_min);

  boil::cart.min_int(&nontrivial_0_min);
  boil::cart.min_int(&nontrivial_1_min);
  boil::cart.min_int(&nontrivial_2_min);

  boil::cart.max_int(&aligned_0_max);
  boil::cart.max_int(&aligned_1_max);
  boil::cart.max_int(&aligned_2_max);

  boil::cart.max_int(&nonaligned_0_max);
  boil::cart.max_int(&nonaligned_1_max);
  boil::cart.max_int(&nonaligned_2_max);

  boil::cart.max_int(&nontrivial_0_max);
  boil::cart.max_int(&nontrivial_1_max);
  boil::cart.max_int(&nontrivial_2_max);

#ifdef DEBUG
  boil::oout<<aligned_0_min<<" "<<aligned_0_max<<" "<<nontrivial_0_min<<" "<<nontrivial_0_max
            <<" "<<nonaligned_0_min<<" "<<nonaligned_0_max<<boil::endl;
  boil::oout<<aligned_1_min<<" "<<aligned_1_max<<" "<<nontrivial_1_min<<" "<<nontrivial_1_max
            <<" "<<nonaligned_1_min<<" "<<nonaligned_1_max<<boil::endl;
  boil::oout<<aligned_2_min<<" "<<aligned_2_max<<" "<<nontrivial_2_min<<" "<<nontrivial_2_max
            <<" "<<nonaligned_2_min<<" "<<nonaligned_2_max<<boil::endl;

  for_aijk(i,j,k) {
    stmp[i][j][k] = tempflag2[i][j][k];
  }
  boil::plot->plot(scp,stmp, "sca-flag", time->current_step());
  //exit(0);
#endif

  /* special cases: no aligned in layer 0 or no non-aligned */
  if(  !boil::realistic(aligned_0_min)
     ||!boil::realistic(nonaligned_0_min)
     ||!boil::realistic(nonaligned_1_min)
     ||!boil::realistic(nonaligned_2_min)
    ) {
    insert_bc_curv_HFnormal(scp,ctangential,cnormal,sig);

  } else {
    /**** step 2: identify, if the natural order
                   aligned-(nontrivial)-nonaligned is upheld ****/
    int Icont_0_start, Icont_1_start, Icont_2_start;
    int Icont_0_end, Icont_1_end, Icont_2_end;

    Icont_0_start = Icont_1_start = Icont_2_start = -boil::unint;
    Icont_0_end = Icont_1_end = Icont_2_end = -boil::unint;

    /* 0 */
    if(!boil::realistic(nontrivial_0_min)||nonaligned_0_min<nontrivial_0_min) {
      Icont_0_start = nonaligned_0_min;
      Icont_0_end   = nonaligned_0_min;
    } else {
      Icont_0_start = nontrivial_0_min;
      Icont_0_end   = nonaligned_0_min;
    }

    /* 1 */
    if(!boil::realistic(nontrivial_1_min)||nonaligned_1_min<nontrivial_1_min) {
      Icont_1_start = nonaligned_1_min;
      Icont_1_end   = nonaligned_1_min;
    } else {
      Icont_1_start = nontrivial_1_min;
      Icont_1_end   = nonaligned_1_min;
    }

    /* 2 */
    if(!boil::realistic(nontrivial_2_min)||nonaligned_2_min<nontrivial_2_min) {
      Icont_2_start = nonaligned_2_min;
      Icont_2_end   = nonaligned_2_min;
    } else {
      Icont_2_start = nontrivial_2_min;
      Icont_2_end   = nonaligned_2_min;
    }

    /* exception: film-like geometry */
    int Nfilm = Icont_0_end-Icont_0_start+1;

#ifdef DEBUG
    boil::oout<<Icont_0_start<<" "<<Icont_0_end<<boil::endl;
    boil::oout<<Icont_1_start<<" "<<Icont_1_end<<boil::endl;
    boil::oout<<Icont_2_start<<" "<<Icont_2_end<<boil::endl;
    boil::oout<<Nfilm<<boil::endl;
#endif

    if(Nfilm>Nfilm_crit) {
      insert_bc_curv_HFnormal(scp,ctangential,cnormal,sig);

    } else {
      /* to find the stencil extent, we need to check the maximum of Icont */
      int Icont_max_end = std::max(Icont_0_end,
                                   std::max(Icont_1_end,Icont_2_end));

      /* now we need to check for possible corrugation failure */
      for_ijk(i,j,k) {
        if(dom->ibody().on(i,j,k)) {
          /* 0 */
          if(dom->ibody().off(i,j,k-1)||(kminw && (k==sk()  ))) {
            int Iglob_0 = scp.domain()->global_I(i)-boil::BW+1;
            if(Iglob_0>=Icont_0_start&&Iglob_0<=Icont_max_end) {
              if(tempflag2[i][j][k]==1) {
                Icont_0_end = boil::unint;
              }
            }
          }

          /* 1 */
          if(dom->ibody().off(i,j,k-2)||(kminw && (k==sk()+1))) {
            int Iglob_1 = scp.domain()->global_I(i)-boil::BW+1;
            if(Iglob_1>=Icont_1_start&&Iglob_1<=Icont_max_end) {
              if(tempflag2[i][j][k]==1) {
                Icont_1_end = boil::unint;
              }
            }
          }

          /* 2 */
          if(dom->ibody().off(i,j,k-3)||(kminw && (k==sk()+2))) {
            int Iglob_2 = scp.domain()->global_I(i)-boil::BW+1;
            if(Iglob_2>=Icont_2_start&&Iglob_2<=Icont_max_end) {
              if(tempflag2[i][j][k]==1) {
                Icont_2_end = boil::unint;
              }
            }
          }

        } /* on */
      } /* ijk */
      boil::cart.max_int(&Icont_0_end);
      boil::cart.max_int(&Icont_1_end);
      boil::cart.max_int(&Icont_2_end);

#ifdef DEBUG
      boil::oout<<Icont_0_end<<" ";
      boil::oout<<Icont_1_end<<" ";
      boil::oout<<Icont_2_end<<" | ";
      boil::oout<<Icont_max_end<<boil::endl;
#endif

      /* if any of them is unreal, geometry is unsuitable for parallel HF */
      if(  !boil::realistic(Icont_0_end)
         ||!boil::realistic(Icont_1_end)  
         ||!boil::realistic(Icont_2_end)  
        ) {
        insert_bc_curv_HFnormal(scp,ctangential,cnormal,sig);

      } else {
        /* contact-line type of boundary condition (yay!) */
        insert_bc_curv_HFparallel(scp,ctangential,cnormal,sig,
                                  Range<int>(1,Icont_max_end));
        if(Icont_max_end<scp.domain()->gi()) {
          insert_bc_curv_HFnormal(scp,ctangential,cnormal,sig,
                                  Range<int>(Icont_max_end+1,scp.domain()->gi()));
        }
      }
    } /* not film */
  } /* not special case */

  return;
  
}
