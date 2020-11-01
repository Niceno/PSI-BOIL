#include "heaviside.h"

/******************************************************************************/
void Heaviside::fs_bnd_nosubgrid(const Scalar & scp, Vector & fs,
                                 const real & tol_wall) {
/***************************************************************************//**
*  \brief Corrects fs at boundaries. No subgrid interfaces are considered.
*         IMPORTANT: does not work when immersed boundaries do not correspond
*                    to cell boundaries!!!
******************************************************************************/

  for( int b=0; b<scp.bc().count(); b++ ) {

    if( scp.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( scp.bc().type(b) == BndType::wall() ) {

      Dir d = scp.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int of(0), ofx(0), ofy(0), ofz(0);
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          of  = +1;
          ofx = +1;
        } else if (d == Dir::imax()) { 
          mcomp = Comp::i();
          of  = 0;
          ofx = -1;
        } else if (d == Dir::jmin()) { 
          mcomp = Comp::j();
          of  = +1;
          ofy = +1;
        } else if (d == Dir::jmax()) { 
          mcomp = Comp::j();
          of  = 0;
          ofy = -1;
        } else if (d == Dir::kmin()) { 
          mcomp = Comp::k();
          of  = +1;
          ofz = +1;
        } else if (d == Dir::kmax()) { 
          mcomp = Comp::k();
          of  = 0;
          ofz = -1;
        } else {
          continue;
        }
       
        for_vijk( scp.bc().at(b), i,j,k ) { 
          /* the fs value is reset */
          if(mcomp==Comp::i()) {
            fs[mcomp][i+of][j][k] = boil::unreal;
          } else if(mcomp==Comp::j()) {
            fs[mcomp][i][j+of][k] = boil::unreal;
          } else {
            fs[mcomp][i][j][k+of] = boil::unreal;
          }
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);
         
    /* erroneous interfaces */
    real scpscp = scp[i][j][k];
    bool errint = (scpscp<tol_wall||scpscp-1.0>-tol_wall);
    if(errint) {
      if(dom->ibody().off(i-1,j,k)) {
        mcomp = Comp::i();
        fs[mcomp][i  ][j][k] = boil::unreal;
      }
      if(dom->ibody().off(i+1,j,k)) {
        mcomp = Comp::i();
        fs[mcomp][i+1][j][k] = boil::unreal;
      }
      if(dom->ibody().off(i,j-1,k)) {
        mcomp = Comp::j();
        fs[mcomp][i][j  ][k] = boil::unreal;
      }
      if(dom->ibody().off(i,j+1,k)) {
        mcomp = Comp::j();
        fs[mcomp][i][j+1][k] = boil::unreal;
      }
      if(dom->ibody().off(i,j,k-1)) {
        mcomp = Comp::k();
        fs[mcomp][i][j][k  ] = boil::unreal;
      }
      if(dom->ibody().off(i,j,k+1)) {
        mcomp = Comp::k();
        fs[mcomp][i][j][k+1] = boil::unreal;
      }
      continue;
    }

    /* west is in solid domain */
    if(dom->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      fs[mcomp][i  ][j][k] = boil::unreal;
    }

    /* east */
    if(dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      fs[mcomp][i+1][j][k] = boil::unreal;
    }

    /* south */
    if(dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      fs[mcomp][i][j  ][k] = boil::unreal;
    }

    /* north */
    if(dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      fs[mcomp][i][j+1][k] = boil::unreal;
    }

    /* bottom */
    if(dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      fs[mcomp][i][j][k  ] = boil::unreal;
    }

    /* top */
    if(dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      fs[mcomp][i][j][k+1] = boil::unreal;
    }
  }

  return;
}
