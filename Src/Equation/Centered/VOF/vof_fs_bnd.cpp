#include "vof.h"

/******************************************************************************/
void VOF::fs_bnd(const Scalar & scp) {
/***************************************************************************//**
*  \brief Corrects fs at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
*         IMPORTANT: does not work when immersed boundaries do not correspond
*                    to cell boundaries!!!
******************************************************************************/

  /* tolerance is necessary because of errors */
  /* e.g. 1.0 approx 0.999 */
  real tolf = 0.0e-2;
  /* tol_wall is defined in header */
  //real tol_wall = 0.5e-2; 
  /* consistent tol_wall with update_at_walls: ie for clr */

  for( int b=0; b<scp.bc().count(); b++ ) {

    if( scp.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( scp.bc().type(b) == BndType::wall() ) {

      Dir d = scp.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int of(0), ofx(0), ofy(0), ofz(0);
#if 0 /* doesn't work without constexpr */
        switch(d) {
          case Dir::imin() : mcomp = Comp::i();
                             of  = +1;
                             ofx = +1;
                             break;  
          case Dir::imax() : mcomp = Comp::i();
                             of  = 0;
                             ofx = -1;
                             break;  
          case Dir::jmin() : mcomp = Comp::j();
                             of  = +1;
                             ofy = +1;
                             break;  
          case Dir::jmax() : mcomp = Comp::j();
                             of  = 0;
                             ofy = -1;
                             break;  
          case Dir::kmin() : mcomp = Comp::k();
                             of  = +1;
                             ofz = +1;
                             break;  
          case Dir::kmax() : mcomp = Comp::k();
                             of  = 0;
                             ofz = -1;
                             break;  
          default : continue;
        }
#else
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
#endif
       
        for_vijk( scp.bc().at(b), i,j,k ) { 
          /* at first, the fs value is reset */
          if(mcomp==Comp::i()) {
            fs[mcomp][i+of][j][k] = boil::unreal;
          } else if(mcomp==Comp::j()) {
            fs[mcomp][i][j+of][k] = boil::unreal;
          } else {
            fs[mcomp][i][j][k+of] = boil::unreal;
          }

          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
          real fsval = fs_val(mcomp,ii,jj,kk);
          bool flagm = (0.0+tolf <= fsval && fsval <= 0.5    );
          bool flagp = (0.5     <= fsval && fsval <= 1.0-tolf);
         
          /* erroneous interfaces */
          real scpscp = scp[ii][jj][kk];
          bool errint = (scpscp<tol_wall||scpscp-1.0>-tol_wall);

          /* interface exists in the dir we are interested in */
          if(((flagm && of)||(flagp && !of))&&!errint) { 
            if(mcomp==Comp::i()) {
              fs[mcomp][i+of][j][k] = scp.xn(ii) + scp.dxc(ii) * fsval;
            } else if(mcomp==Comp::j()) {
              fs[mcomp][i][j+of][k] = scp.yn(jj) + scp.dyc(jj) * fsval;   
            } else {
              fs[mcomp][i][j][k+of] = scp.zn(kk) + scp.dzc(kk) * fsval;
            }
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
    if (dom->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tolf <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i  ][j][k] = scp.xn(i) + scp.dxc(i) * fsval;
      else
        fs[mcomp][i  ][j][k] = boil::unreal;
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tolf);
      if(flagp)
        fs[mcomp][i+1][j][k] = scp.xn(i) + scp.dxc(i) * fsval;
      else
        fs[mcomp][i+1][j][k] = boil::unreal;
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tolf <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j  ][k] = scp.yn(j) + scp.dyc(j) * fsval;
      else
        fs[mcomp][i][j  ][k] = boil::unreal;
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tolf);
      if(flagp)
        fs[mcomp][i][j+1][k] = scp.yn(j) + scp.dyc(j) * fsval;
      else
        fs[mcomp][i][j+1][k] = boil::unreal;
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tolf <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j][k  ] = scp.zn(k) + scp.dzc(k) * fsval;
      else
        fs[mcomp][i][j][k  ] = boil::unreal;
      
      //boil::oout<<"VOF-fs_bnd "<<i<<" "<<j<<" "<<k<<" | "<<scp.zn(k)<<" "<<fsval<<" "<<scp.dzc(k)<<" "<<fs[mcomp][i][j][k  ]<<boil::endl;
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tolf);
      if(flagp)
        fs[mcomp][i][j][k+1] = scp.zn(k) + scp.dzc(k) * fsval;
      else
        fs[mcomp][i][j][k+1] = boil::unreal;
    }
  }

  return;
}
