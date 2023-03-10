#include "topology.h"

/******************************************************************************/
void Topology::fs_bnd_geometric(const Scalar & scp, Vector & fs, 
                                const real & tol_wall) {
/***************************************************************************//**
*  \brief Corrects fs at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
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
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
         
          /* erroneous interfaces */
          real scpscp = scp[ii][jj][kk];
          bool errint = (scpscp<tol_wall||scpscp-1.0>-tol_wall);

          /* the fs value is reset if erroneous or if it is inside wall */
          if(mcomp==Comp::i()) {
            real & fsval = fs[mcomp][i+of][j][k];
            real refpos = scp.xn(i+of);
            bool test = (d == Dir::imin()) ? (fsval <= refpos) : (fsval >= refpos);
            if(errint||test)
              fsval = boil::unreal;
          } else if(mcomp==Comp::j()) {
            real & fsval = fs[mcomp][i][j+of][k];
            real refpos = scp.yn(j+of);
            bool test = (d == Dir::jmin()) ? (fsval <= refpos) : (fsval >= refpos);
            if(errint||test)
              fsval = boil::unreal;
          } else {
            real & fsval = fs[mcomp][i][j][k+of];
            real refpos = scp.zn(k+of);
            bool test = (d == Dir::kmin()) ? (fsval <= refpos) : (fsval >= refpos);
            if(errint||test)
              fsval = boil::unreal;
          }
        }

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<domain()->ibody().nccells(); cc++) {
    int i,j,k;
    Comp mcomp;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    domain()->ibody().ijk(cc,&i,&j,&k);
         
    /* erroneous interfaces */
    real scpscp = scp[i][j][k];
    bool errint = (scpscp<tol_wall||scpscp-1.0>-tol_wall);
    if(errint) {
      if(domain()->ibody().off(i-1,j,k)) {
        mcomp = Comp::i();
        fs[mcomp][i  ][j][k] = boil::unreal;
      }
      if(domain()->ibody().off(i+1,j,k)) {
        mcomp = Comp::i();
        fs[mcomp][i+1][j][k] = boil::unreal;
      }
      if(domain()->ibody().off(i,j-1,k)) {
        mcomp = Comp::j();
        fs[mcomp][i][j  ][k] = boil::unreal;
      }
      if(domain()->ibody().off(i,j+1,k)) {
        mcomp = Comp::j();
        fs[mcomp][i][j+1][k] = boil::unreal;
      }
      if(domain()->ibody().off(i,j,k-1)) {
        mcomp = Comp::k();
        fs[mcomp][i][j][k  ] = boil::unreal;
      }
      if(domain()->ibody().off(i,j,k+1)) {
        mcomp = Comp::k();
        fs[mcomp][i][j][k+1] = boil::unreal;
      }
      continue;
    }

    /* west is in solid domain */
    if(domain()->ibody().off(i-1,j,k)) {
      mcomp = Comp::i();
      real & fsval = fs[mcomp][i  ][j][k];
      real refpos = scp.xn(i);
      bool test = (fsval <= refpos);
      if(test)
        fsval = boil::unreal;
    }

    /* east */
    if(domain()->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      real & fsval = fs[mcomp][i+1][j][k];
      real refpos = scp.xn(i+1);
      bool test = (fsval >= refpos);
      if(test)
        fsval = boil::unreal;
    }

    /* south */
    if(domain()->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      real & fsval = fs[mcomp][i][j  ][k];
      real refpos = scp.yn(j);
      bool test = (fsval <= refpos);
      if(test)
        fsval = boil::unreal;
    }

    /* north */
    if(domain()->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      real & fsval = fs[mcomp][i][j+1][k];
      real refpos = scp.yn(j+1);
      bool test = (fsval >= refpos);
      if(test)
        fsval = boil::unreal;
    }

    /* bottom */
    if(domain()->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      real & fsval = fs[mcomp][i][j][k  ];
      real refpos = scp.zn(k);
      bool test = (fsval <= refpos);
      if(test)
        fsval = boil::unreal;
    }

    /* top */
    if(domain()->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      real & fsval = fs[mcomp][i][j][k+1];
      real refpos = scp.zn(k+1);
      bool test = (fsval >= refpos);
      if(test)
        fsval = boil::unreal;
    }
  }

  return;
}
