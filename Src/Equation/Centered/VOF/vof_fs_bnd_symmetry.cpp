#include "vof.h"

/******************************************************************************/
void VOF::fs_bnd_symmetry(const Scalar & scp, Vector & fs, 
                          const real & tol_wall) {
/***************************************************************************//**
*  \brief fs at symmetry plane.
******************************************************************************/

  for( int b=0; b<scp.bc().count(); b++ ) {

    if( scp.bc().type_decomp(b) ) continue;

    /* special treatment at walls */
    if( scp.bc().type(b) == BndType::symmetry() ) {

      Dir d = scp.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int ofx(0), ofy(0), ofz(0);
        int oofx(0), oofy(0), oofz(0);
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          oofx= +1;
          ofx = +1;
        } else if (d == Dir::imax()) { 
          mcomp = Comp::i();
          oofx= 0;
          ofx = -1;
        } else if (d == Dir::jmin()) { 
          mcomp = Comp::j();
          oofy= +1;
          ofy = +1;
        } else if (d == Dir::jmax()) { 
          mcomp = Comp::j();
          oofy= 0;
          ofy = -1;
        } else if (d == Dir::kmin()) { 
          mcomp = Comp::k();
          oofz= +1;
          ofz = +1;
        } else if (d == Dir::kmax()) { 
          mcomp = Comp::k();
          oofz= 0;
          ofz = -1;
        } else {
          continue;
        }
       
        for_vijk( scp.bc().at(b), i,j,k ) { 

          /* step one: mirror existing values */
          int iref(i+oofx),jref(j+oofy),kref(k+oofz);

          for(int bf(1); bf <= boil::BW; ++bf) {
            int ii_in = iref+bf*ofx;
            int jj_in = jref+bf*ofy;
            int kk_in = kref+bf*ofz;
            int ii_out = iref-bf*ofx;
            int jj_out = jref-bf*ofy;
            int kk_out = kref-bf*ofz;
            
            real refpos = (mcomp==Comp::i())
                            ? scp.xn(iref)    /* if i() */
                        : (mcomp==Comp::j())
                            ? scp.yn(jref)    /* if j() */
                            : scp.zn(kref);   /* if k() */

            real fspos = fs[mcomp][ii_in][jj_in][kk_in];
            if(boil::realistic(fspos))
              fs[mcomp][ii_out][jj_out][kk_out] = 2.*refpos - fspos;
          }

          /* step two: set values in near-bnd cells */
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;
         
          /* erroneous interfaces */
          real scpscp = scp[ii][jj][kk];
          bool errint = (scpscp<tol_wall||scpscp-1.0>-tol_wall);

          if(errint)
            continue;

          /* the fs value is reset if erroneous or if it is inside wall */
          if(mcomp==Comp::i()) {
            real fsval = fs_val(mcomp,ii,jj,kk);
            real refpos1 = scp.xn(i+oofx);
            real refpos2 = scp.xc(ii);
            real refpos3 = scp.xn(ii);
            real fspos = refpos3 + scp.dxc(ii)*fsval;
            bool test = (d == Dir::imin()) 
                      ? ((fspos <= refpos1)||(fspos >= refpos2))
                      : ((fspos >= refpos1)||(fspos <= refpos2));
            if(!test)
              fs[mcomp][i+oofx][j][k] = fspos;

          } else if(mcomp==Comp::j()) {
            real fsval = fs_val(mcomp,ii,jj,kk);
            real refpos1 = scp.yn(j+oofy);
            real refpos2 = scp.yc(jj);
            real refpos3 = scp.yn(jj);
            real fspos = refpos3 + scp.dyc(jj)*fsval;
            bool test = (d == Dir::jmin())
                      ? ((fspos <= refpos1)||(fspos >= refpos2))
                      : ((fspos >= refpos1)||(fspos <= refpos2));
            if(!test)
              fs[mcomp][i][j+oofy][k] = fspos;

          } else {
            real fsval = fs_val(mcomp,ii,jj,kk);
            real refpos1 = scp.zn(k+oofz);
            real refpos2 = scp.zc(kk);
            real refpos3 = scp.zn(kk);
            real fspos = refpos3 + scp.dzc(kk)*fsval;
            bool test = (d == Dir::kmin())
                      ? ((fspos <= refpos1)||(fspos >= refpos2))
                      : ((fspos >= refpos1)||(fspos <= refpos2));
            if(!test)
              fs[mcomp][i][j][k+oofz] = fspos;
          }
        } /* vijk */

      } /* dir not undefined */
    } /* is symmetry? */
  } /* loop over bcs */

}
