#include "finescalar.h"

/******************************************************************************/
void FineScalar::bdcond_edge() {
/***************************************************************************//**
*  \brief Correct the marker function at edges at boundaries
*******************************************************************************/

  /* simplified way, assuming fs was already evaluated */

  for( int bb=0; bb<phi->bc().count(); bb++ ) {

    if( phi->bc().type_decomp(bb) ) continue;

    /* special treatment at walls */
    if( phi->bc().type(bb) == BndType::wall() ) {

      Dir d = phi->bc().direction(bb);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int ofxv(0), ofyv(0), ofzv(0), ofxs(0), ofys(0), ofzs(0),dir;
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          dir = w();
          ofxv = +1;
          ofxs = +1;
        } else if (d == Dir::imax()) {
          mcomp = Comp::i();
          dir = e();
          ofxv = 0;
          ofxs = -1;
        } else if (d == Dir::jmin()) {
          mcomp = Comp::j();
          dir = s();
          ofyv = +1;
          ofys = +1;
        } else if (d == Dir::jmax()) {
          mcomp = Comp::j();
          dir = n();
          ofyv = 0;
          ofys = -1;
        } else if (d == Dir::kmin()) {
          mcomp = Comp::k();
          dir = b();
          ofzv = +1;
          ofzs = +1;
        } else if (d == Dir::kmax()) {
          mcomp = Comp::k();
          dir = t();
          ofzv = 0;
          ofzs = -1;
        } else {
          continue;
        }

        for_vijk( phi->bc().at(bb), i,j,k ) {

          int ii = i+ofxs;
          int jj = j+ofys;
          int kk = k+ofzs;

          real fsval = (*fs)[mcomp][i+ofxv][j+ofyv][k+ofzv];
          real phiref = (*phi)[ii][jj][kk];

          if(fsval>boil::zetta) {
            value(ii,jj,kk,dir) = phiref>phisurf;
          } else {
            value(ii,jj,kk,dir) = phiref<phisurf;
          }
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*--------------+
  | immersed body |
  +--------------*/
  const Domain * dom = phi->domain();
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    int dir;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in solid domain */
    if (dom->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      dir = w();
      real fsval = (*fs)[mcomp][i  ][j][k];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      dir = e();
      real fsval = (*fs)[mcomp][i+1][j][k];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      dir = s();
      real fsval = (*fs)[mcomp][i][j  ][k];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      dir = n();
      real fsval = (*fs)[mcomp][i][j+1][k];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      dir = b();
      real fsval = (*fs)[mcomp][i][j][k  ];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      dir = t();
      real fsval = (*fs)[mcomp][i][j][k+1];
      real phiref = (*phi)[i][j][k];

      if(fsval>boil::zetta) {
        value(i,j,k,dir) = phiref>phisurf;
      } else {
        value(i,j,k,dir) = phiref<phisurf;
      }
    }
  }

}

