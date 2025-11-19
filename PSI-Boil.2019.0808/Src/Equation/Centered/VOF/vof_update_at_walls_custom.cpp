#include "vof.h"

/******************************************************************************/
void VOF::update_at_walls_custom(Scalar & scp) {
/***************************************************************************//**
*  \brief Prepares volume fraction for marching cube at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
******************************************************************************/

  /* tolerance is necessary because of errors */
  /* e.g. 1.0 approx 0.999 */
  /* now defined in the header file */
  //real tol_wall = 0.5e-2;

  /*-------------+
  | single walls |
  +-------------*/
  for( int b=0; b<scp.bc().count(); b++ ) {

    if( scp.bc().type_decomp(b) ) continue;

    if( scp.bc().type(b) == BndType::wall() ) {

      Sign sig;
      Comp mcomp;
      Dir d = scp.bc().direction(b);
      if(d != Dir::undefined()) {
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          sig = Sign::pos();
        } else if (d == Dir::imax()) { 
          mcomp = Comp::i();
          sig = Sign::neg();
        } else if (d == Dir::jmin()) { 
          mcomp = Comp::j();
          sig = Sign::pos();
        } else if (d == Dir::jmax()) { 
          mcomp = Comp::j();
          sig = Sign::neg();
        } else if (d == Dir::kmin()) { 
          mcomp = Comp::k();
          sig = Sign::pos();
        } else if (d == Dir::kmax()) { 
          mcomp = Comp::k();
          sig = Sign::neg();
        } else {
          continue;
        }
       
        for_vijk( scp.bc().at(b), i,j,k ) { 
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*--------------+
  | immersed body |
  +--------------*/
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) {

        Sign sig;
        Comp mcomp;

        /* west is in fluid domain */
        if(dom->ibody().on(i-1,j,k)) {
          mcomp = Comp::i();
          sig = Sign::neg();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

        /* east is in fluid domain */
        if(dom->ibody().on(i+1,j,k)) {
          mcomp = Comp::i();
          sig = Sign::pos();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

        /* south is in fluid domain */
        if(dom->ibody().on(i,j-1,k)) {
          mcomp = Comp::j();
          sig = Sign::neg();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

        /* north is in fluid domain */
        if(dom->ibody().on(i,j+1,k)) {
          mcomp = Comp::j();
          sig = Sign::pos();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

        /* bottom is in fluid domain */
        if(dom->ibody().on(i,j,k-1)) {
          mcomp = Comp::k();
          sig = Sign::neg();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

        /* top is in fluid domain */
        if(dom->ibody().on(i,j,k+1)) {
          mcomp = Comp::k();
          sig = Sign::pos();
          scp[i][j][k] = update_at_walls_func(sig,mcomp,i,j,k);
        }

      } /* center off */
    } /* ijk */
  } /* ibody exists */

  scp.bnd_update_nowall();
  scp.exchange_all();
  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();
 
  return;
}
