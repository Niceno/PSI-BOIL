#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::update_at_walls(Scalar & sca) {
/***************************************************************************//**
*  \brief Prepares color function for marching cube at boundaries.
*         scalar_exchange(_all) should take account of periodic condition.
*         IMPORTANT: does not work when immersed boundaries do not correspond
*                    to cell boundaries!!!
*         MOREOVER:  since the normal vector is not properly corrected at
*                    boundaries consisting partially of walls and partially
*                    of immersed bodies, this special case was not considered.
*                    if such a thing was to be implemented, cases like
*                    ib-wall edge, wall-ib edge, wall-ib-wall corner,
*                    wall-ib-ib corner etc would have to be included!!!
******************************************************************************/

  /* tolerance is necessary because of errors */
  /* e.g. 1.0 approx 0.999 */
  /* now defined in the header file */
  //real tol_wall = 0.5e-2;

  /*-------------+
  | single walls |
  +-------------*/
  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    if( sca.bc().type(b) == BndType::wall() ) {

      Dir d = sca.bc().direction(b);
      if(d != Dir::undefined()) {
        int ofx(0),ofy(0),ofz(0);
        real rat;
        
        if       (d == Dir::imin()) {
          ofx = +1;
          rat = -sca.dxe(si())/sca.dxc(si());
        } else if(d == Dir::imax()) { 
          ofx = -1;
          rat = -sca.dxw(ei())/sca.dxc(ei());
        } else if(d == Dir::jmin()) { 
          ofy = +1;
          rat = -sca.dyn(sj())/sca.dyc(sj());
        } else if(d == Dir::jmax()) { 
          ofy = -1;
          rat = -sca.dys(ej())/sca.dyc(ej());
        } else if(d == Dir::kmin()) { 
          ofz = +1;
          rat = -sca.dzt(sk())/sca.dzc(sk());
        } else if(d == Dir::kmax()) { 
          ofz = -1;
          rat = -sca.dzb(ek())/sca.dzc(ek());
        } else {
          continue;
        }

        for_vijk( sca.bc().at(b), i,j,k ) { 
          sca[i][j][k] = extrapolate_c(sca,i,j,k,ofx,ofy,ofz,rat);
        } /* for loop */
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

#if 0 /* underdevelopment */
  /*---------------+
  | edges of walls |
  +---------------*/
  
  /* line i-min & j-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=si()-1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = +1;
    xpos = -0.5;
    ypos = -0.5;
    for_k(k) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & j-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=si()-1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = -1;
    xpos = -0.5;
    ypos = 1.5;
    for_k(k) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & j-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = +1;
    xpos = 1.5;
    ypos = -0.5;
    for_k(k) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & j-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = -1;
    xpos = 1.5;
    ypos = 1.5;
    for_k(k) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = +1;
    xpos = -0.5;
    zpos = -0.5;
    for_j(j) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = -1;
    xpos = -0.5;
    zpos = 1.5;
    for_j(j) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofz = +1;
    xpos = 1.5;
    zpos = -0.5;
    for_j(j) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_j(j) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-min & k-min */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=sj()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = +1;
    ypos = -0.5;
    zpos = -0.5;
    for_i(i) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-min & k-max */
  if(sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=sj()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = -1;
    ypos = -0.5;
    zpos = 1.5;
    for_i(i) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-max & k-min */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = +1;
    ypos = 1.5;
    zpos = -0.5;
    for_i(i) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-max & k-max */
  if(sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_i(i) {
      sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /*-----------------+
  | corners of walls |
  +-----------------*/

  /* corner i-min & j-min & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si()-1;
    int j=sj()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = +1;
    ofz = +1;
    xpos = -0.5;
    ypos = -0.5;
    zpos = -0.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-min & j-min & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si()-1;
    int j=sj()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = +1;
    ofz = -1;
    xpos = -0.5;
    ypos = -0.5;
    zpos = 1.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-min & j-max & k-min */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si()-1;
    int j=ej()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = -1;
    ofz = +1;
    xpos = -0.5;
    ypos = 1.5;
    zpos = -0.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-min & j-max & k-max */
  if(sca.bc().type_here(Dir::imin(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si()-1;
    int j=ej()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = -1;
    ofz = -1;
    xpos = -0.5;
    ypos = 1.5;
    zpos = 1.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-max & j-min & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=sj()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = +1;
    ofz = +1;
    xpos = 1.5;
    ypos = -0.5;
    zpos = -0.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-min & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmin(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=sj()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = +1;
    ofz = -1;
    xpos = 1.5;
    ypos = -0.5;
    zpos = 1.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-max & k-min */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = -1;
    ofz = +1;
    xpos = 1.5;
    ypos = 1.5;
    zpos = -0.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-max & k-max */
  if(sca.bc().type_here(Dir::imax(), BndType::wall()) &&
     sca.bc().type_here(Dir::jmax(), BndType::wall()) &&
     sca.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = -1;
    ofz = -1;
    xpos = 1.5;
    ypos = 1.5;
    zpos = 1.5;
    sca[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }
#endif

  /*--------------+
  | immersed body |
  +--------------*/
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in solid domain */
    if(dom->ibody().off(i-1,j,k)) {
      int ofx(0), ofy(0), ofz(0);
      ofx = +1;
      real rat = -sca.dxe(i)/sca.dxc(i);
      sca[i-1][j][k] = extrapolate_c(sca,i-1,j,k,ofx,ofy,ofz,rat);
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      int ofx(0), ofy(0), ofz(0);
      ofx = -1;
      real rat = -sca.dxw(i)/sca.dxc(i);
      sca[i+1][j][k] = extrapolate_c(sca,i+1,j,k,ofx,ofy,ofz,rat);
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      int ofx(0), ofy(0), ofz(0);
      ofy = +1;
      real rat = -sca.dyn(j)/sca.dyc(j);
      sca[i][j-1][k] = extrapolate_c(sca,i,j-1,k,ofx,ofy,ofz,rat);
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      int ofx(0), ofy(0), ofz(0);
      ofy = -1;
      real rat = -sca.dys(j)/sca.dyc(j);
      sca[i][j+1][k] = extrapolate_c(sca,i,j+1,k,ofx,ofy,ofz,rat);
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      int ofx(0), ofy(0), ofz(0);
      ofz = +1;
      real rat = -sca.dzt(k)/sca.dzc(k);
      sca[i][j][k-1] = extrapolate_c(sca,i,j,k-1,ofx,ofy,ofz,rat);
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      int ofx(0), ofy(0), ofz(0);
      ofy = -1;
      real rat = -sca.dzb(k)/sca.dzc(k);
      sca[i][j][k+1] = extrapolate_c(sca,i,j,k+1,ofx,ofy,ofz,rat);
    }
  }

  sca.exchange_all();
 
  return;
}

/*-------------------+
| ancillary function |
+-------------------*/
real CIPCSL2::extrapolate_c(const Scalar & sca,
                            const int i, const int j, const int k,
                            const int ofx, const int ofy, const int ofz,
                            const real rat) {

  int ii1 = i+ofx;
  int jj1 = j+ofy;
  int kk1 = k+ofz;
 
  int ii2 = i+2*ofx;
  int jj2 = j+2*ofy;
  int kk2 = k+2*ofz;
   

  /* the picos are to avoid singularity for atanh */
  real sca1 = std::max(boil::pico,std::min(1.-boil::pico,sca[ii1][jj1][kk1]));
  real sca2 = std::max(boil::pico,std::min(1.-boil::pico,sca[ii2][jj2][kk2]));

  /* erroneous interfaces */
  if(sca1<tol_wall||sca1-1.0>-tol_wall) {
    return real(sca1>phisurf);
  } 

  real val1 = atanh(2.*sca1-1.); 
  real val2 = atanh(2.*sca2-1.); 
  real arg = val1 + rat*(val2-val1);
 
  return 0.5*(1.+tanh(arg));

}
