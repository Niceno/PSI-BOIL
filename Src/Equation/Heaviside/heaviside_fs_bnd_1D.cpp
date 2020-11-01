#include "heaviside.h"

/******************************************************************************/
void Heaviside::fs_bnd_1D(const Scalar & scp, Vector & fs, 
                          const real & tol_wall, const Sign & sig) {
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
          
          bool test = (sig == Sign::pos()) ? (scpscp >= clrsurf) : (scpscp <= clrsurf);

          /* the fs value is reset if erroneous or centre of cell in
             liquid/gas (based on sign) */
          /* otherwise, it is corrected based on SLIC approximation */
          if(mcomp==Comp::i()) {
            real & fsval = fs[mcomp][i+of][j][k];
            real refpos = scp.xn(i+of);
            if(errint||test) {
              fsval = boil::unreal;
            } else {
              real sval = std::min(1.0,std::max(0.0,scpscp));
              if(sig==Sign::neg())
                sval = 1.0-sval;
              if(d==Dir::imax())
                sval = 1.0-sval;
              fsval = scp.xn(ii) + sval*scp.dxc(ii);
            }
          } else if(mcomp==Comp::j()) {
            real & fsval = fs[mcomp][i][j+of][k];
            real refpos = scp.yn(j+of);
            if(errint||test) {
              fsval = boil::unreal;
            } else {
              real sval = std::min(1.0,std::max(0.0,scpscp));
              if(sig==Sign::neg())
                sval = 1.0-sval;
              if(d==Dir::jmax())
                sval = 1.0-sval;
              fsval = scp.yn(jj) + sval*scp.dyc(jj);
            }
          } else {
            real & fsval = fs[mcomp][i][j][k+of];
            real refpos = scp.zn(k+of);
            if(errint||test) {
              fsval = boil::unreal;
            } else {
              real sval = std::min(1.0,std::max(0.0,scpscp));
              if(sig==Sign::neg())
                sval = 1.0-sval;
              if(d==Dir::kmax())
                sval = 1.0-sval;
              fsval = scp.zn(kk) + sval*scp.dzc(kk);
            }
          }
        } /* for vijk */

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
    bool test = (sig == Sign::pos()) ? (scpscp >= clrsurf) : (scpscp <= clrsurf);
    if(errint||test) {
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
      real & fsval = fs[mcomp][i  ][j][k];
      real sval = std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.xn(i) + sval*scp.dxc(i);
    }

    /* east */
    if(dom->ibody().off(i+1,j,k)) {
      mcomp = Comp::i();
      real & fsval = fs[mcomp][i+1][j][k];
      real sval = 1.0 - std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.xn(i) + sval*scp.dxc(i);
    }

    /* south */
    if(dom->ibody().off(i,j-1,k)) {
      mcomp = Comp::j();
      real & fsval = fs[mcomp][i][j  ][k];
      real sval = std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.yn(j) + sval*scp.dyc(j);
    }

    /* north */
    if(dom->ibody().off(i,j+1,k)) {
      mcomp = Comp::j();
      real & fsval = fs[mcomp][i][j+1][k];
      real sval = 1.0 - std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.yn(j) + sval*scp.dyc(j);
    }

    /* bottom */
    if(dom->ibody().off(i,j,k-1)) {
      mcomp = Comp::k();
      real & fsval = fs[mcomp][i][j][k  ];
      real sval = std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.zn(k) + sval*scp.dzc(k);
    }

    /* top */
    if(dom->ibody().off(i,j,k+1)) {
      mcomp = Comp::k();
      real & fsval = fs[mcomp][i][j][k+1];
      real sval = 1.0 - std::min(1.0,std::max(0.0,scpscp));
      if(sig==Sign::neg())
        sval = 1.0-sval;
      fsval = scp.zn(k) + sval*scp.dzc(k);
    }
  }

  return;
}
