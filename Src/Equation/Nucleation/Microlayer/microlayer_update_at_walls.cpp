#include "microlayer.h"

/******************************************************************************/
void Microlayer::update_at_walls(Scalar & clr, Vector & fs) { 
/***************************************************************************//**
*  \brief Corrects clr & fs due to subgrid microlayer existence.
******************************************************************************/

  for( int b=0; b<dmicro.bc().count(); b++ ) {

    if( dmicro.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( dmicro.bc().type(b) == BndType::wall() ) {

      Dir d = dmicro.bc().direction(b);
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
       
        for_vijk( dmicro.bc().at(b), i,j,k ) { 
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
         
          real scpscp = dmicro[ii][jj][kk];
          if ( scpscp <= dmicro_min*(1.+boil::pico) // depleted
            || !boil::realistic(scpscp)) {          // full vapor
          } else {
            clr[i][j][k] = 1.0;
            //(*vf)[ii][jj][kk] = scpscp/(2.*cht->distance_face(Sign::neg(),mcomp,ii,jj,kk));
            //(*cht->topo->clr)[ii][jj][kk] = scpscp/(2.*cht->distance_face(Sign::neg(),mcomp,ii,jj,kk));
            if(mcomp==Comp::i()) {
              real & fsval = fs[mcomp][i+of][j][k];
              real refpos = dmicro.xn(i+of);
              real sval = scpscp/dmicro.dxc(ii); //1.1*tol_wall;
              if(d==Dir::imax())
                sval = 1.0-sval;
              fsval = dmicro.xn(ii) + sval*dmicro.dxc(ii);
            } else if(mcomp==Comp::j()) {
              real & fsval = fs[mcomp][i][j+of][k];
              real refpos = dmicro.yn(j+of);
              real sval = scpscp/dmicro.dyc(jj); //1.1*tol_wall;
              if(d==Dir::jmax())
                sval = 1.0-sval;
              fsval = dmicro.yn(jj) + sval*dmicro.dyc(jj);
            } else {
              real & fsval = fs[mcomp][i][j][k+of];
              real refpos = dmicro.zn(k+of);
              real sval = scpscp/dmicro.dzc(kk); //1.1*tol_wall;
              if(d==Dir::kmax())
                sval = 1.0-sval;
              fsval = dmicro.zn(kk) + sval*dmicro.dzc(kk);
            }
          } /* microlayer is real */
        } /* for vijk */

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<cht->topo->domain()->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    cht->topo->domain()->ibody().ijk(cc,&i,&j,&k);
         
    real scpscp = dmicro[i][j][k];
    if ( scpscp <= dmicro_min*(1.+boil::pico) // depleted
      || !boil::realistic(scpscp)) {          // full vapor
    } else {
      /* west is in solid domain */
      if(cht->topo->domain()->ibody().off(i-1,j,k)) {
        clr[i-1][j][k] = 1.0;
        mcomp = Comp::i();
        real & fsval = fs[mcomp][i  ][j][k];
        real sval = scpscp/dmicro.dxc(i); //1.1*tol_wall;
        fsval = dmicro.xn(i) + sval*dmicro.dxc(i);
      }

      /* east */
      if(cht->topo->domain()->ibody().off(i+1,j,k)) {
        clr[i+1][j][k] = 1.0;
        mcomp = Comp::i();
        real & fsval = fs[mcomp][i+1][j][k];
        real sval = 1.0 - scpscp/dmicro.dxc(i); //1.1*tol_wall;
        fsval = dmicro.xn(i) + sval*dmicro.dxc(i);
      }

      /* south */
      if(cht->topo->domain()->ibody().off(i,j-1,k)) {
        clr[i][j-1][k] = 1.0;
        mcomp = Comp::j();
        real & fsval = fs[mcomp][i][j  ][k];
        real sval = scpscp/dmicro.dyc(j); //1.1*tol_wall;
        fsval = dmicro.yn(j) + sval*dmicro.dyc(j);
      }

      /* north */
      if(cht->topo->domain()->ibody().off(i,j+1,k)) {
        clr[i][j+1][k] = 1.0;
        mcomp = Comp::j();
        real & fsval = fs[mcomp][i][j+1][k];
        real sval = 1.0 - scpscp/dmicro.dyc(j); //1.1*tol_wall;
        fsval = dmicro.yn(j) + sval*dmicro.dyc(j);
      }

      /* bottom */
      if(cht->topo->domain()->ibody().off(i,j,k-1)) {
        clr[i][j][k-1] = 1.0;
        mcomp = Comp::k();
        real & fsval = fs[mcomp][i][j][k  ];
        real sval = scpscp/dmicro.dzc(k); //1.1*tol_wall;
        fsval = dmicro.zn(k) + sval*dmicro.dzc(k);
      }

      /* top */
      if(cht->topo->domain()->ibody().off(i,j,k+1)) {
        clr[i][j][k+1] = 1.0;
        mcomp = Comp::k();
        real & fsval = fs[mcomp][i][j][k+1];
        real sval = 1.0 - scpscp/dmicro.dzc(k); //1.1*tol_wall;
        fsval = dmicro.zn(k) + sval*dmicro.dzc(k);
      }

      //(*vf)[i][j][k] = scpscp/(2.*cht->distance_face(Sign::neg(),mcomp,i,j,k));
      //(*cht->topo->clr)[i][j][k] = scpscp/(2.*cht->distance_face(Sign::neg(),mcomp,i,j,k));

    } /* microlayer real */
  } /* cc */

  clr.exchange_all();

  return;
}
