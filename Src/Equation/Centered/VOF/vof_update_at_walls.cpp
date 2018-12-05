#include "vof.h"

/******************************************************************************/
void VOF::fs_bnd() {
/***************************************************************************//**
*  \brief Prepares volume fraction for marching cube at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
*         IMPORTANT: does not work when immersed boundaries do not correspond
*                    to cell boundaries!!!
******************************************************************************/

  /* tolerance is necessary because of errors */
  /* e.g. 1.0 approx 0.999 */
  real tol = 0.5e-2;

  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( phi.bc().type(b) == BndType::wall() ) {

      Dir d = phi.bc().direction(b);
      if(d != Dir::undefined()) {
        int ofx(0), ofy(0), ofz(0);
        real xpos(0.5), ypos(0.5), zpos(0.5);
        if (d == Dir::imin()) {
          ofx = +1;
          xpos = -0.5; /* below */
        } else if (d == Dir::imax()) { 
          ofx = -1;
          xpos = 1.5;  /* above */
        } else if (d == Dir::jmin()) { 
          ofy = +1;
          ypos = -0.5; /* below */
        } else if (d == Dir::jmax()) { 
          ofy = -1;
          ypos = 1.5;  /* above */
        } else if (d == Dir::kmin()) { 
          ofz = +1;
          zpos = -0.5; /* below */
        } else if (d == Dir::kmax()) { 
          ofz = -1;
          zpos = 1.5;  /* above */
        } else {
          continue;
        }
       
        for_vijk( phi.bc().at(b), i,j,k ) { 
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
   
          real phival = phi[ii][jj][kk];

          /* erroneous interfaces */
          if(phival<tol||phival-1.0>-tol) {
            phi[i][j][k] = real(phival>phisurf);
            continue;
          } 

          /* unnormalized alpha value */
          real alphaval = nalpha[ii][jj][kk];
          
          /* degenerate case I */
          if(!realistic(alphaval)) {
            phi[i][j][k] = real(phival>phisurf);
            continue;
          } 

          /* calculate vn1, vn2, vn3: normal vector at cell center */
          /* n points to the liquid */
          real vn1 = -nx[ii][jj][kk];
          real vn2 = -ny[ii][jj][kk];
          real vn3 = -nz[ii][jj][kk];

          /* surface normal is not affected by translation */
          nx[i][j][k] = -vn1;
          ny[i][j][k] = -vn2;
          nz[i][j][k] = -vn3;
           
          /* normalized normal */
          real vm1 = fabs(vn1);
          real vm2 = fabs(vn2);
          real vm3 = fabs(vn3);

          real denom = vm1+vm2+vm3;

          /* degenerate case II */
          if(denom<boil::pico) {
            phi[i][j][k] = real(phival>phisurf);
            continue;
          }
         
          /* mirror boundary cell to the normalized space */
          if(vn1<0) {
            xpos = 1.0-xpos;
          }
          if(vn2<0) {
            ypos = 1.0-ypos;
          }
          if(vn3<0) {
            zpos = 1.0-zpos;
          }

          /* now, we need to translate the coordinate system */
          /* so that x'(center of bnd cell) = (0.5,0.5,0.5)  */
          /* x' = x - x(center of bnd cell) + (0.5,0.5,0.5)  */
          /* this affects alpha value: */
          /* m dot x' = alpha + m dot [(0.5,0.5,0.5) - x(cbc)] = alpha' */

          alphaval += vm1*(0.5-xpos) + vm2*(0.5-ypos) + vm3*(0.5-zpos);

          /* normalized alpha value */
          alphaval /= denom; 

          /* volume fraction */
          phi[i][j][k] = calc_v(alphaval,vm1,vm2,vm3);
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

#if 0
  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in solid domain */
    if (dom->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i  ][j][k] = phi.xn(i) + phi.dxc(i) * fsval;
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i+1][j][k] = phi.xn(i) + phi.dxc(i) * fsval;
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j  ][k] = phi.yn(j) + phi.dyc(j) * fsval;
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i][j+1][k] = phi.yn(j) + phi.dyc(j) * fsval;
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagm = (0.0+tol <= fsval && fsval <= 0.5    );
      if(flagm)
        fs[mcomp][i][j][k  ] = phi.zn(k) + phi.dzc(k) * fsval;
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      real fsval = fs_val(mcomp,i,j,k);
      bool flagp = (0.5     <= fsval && fsval <= 1.0-tol);
      if(flagp)
        fs[mcomp][i][j][k+1] = phi.zn(k) + phi.dzc(k) * fsval;
    }
  }
#endif
}
