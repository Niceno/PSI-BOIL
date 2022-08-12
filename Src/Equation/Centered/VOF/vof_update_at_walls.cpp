#include "vof.h"

/******************************************************************************/
void VOF::update_at_walls(Scalar & scp) {
/***************************************************************************//**
*  \brief Prepares volume fraction for marching cube at boundaries.
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
  for( int b=0; b<scp.bc().count(); b++ ) {

    if( scp.bc().type_decomp(b) ) continue;

    if( scp.bc().type(b) == BndType::wall() ) {

      Dir d = scp.bc().direction(b);
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
       
        for_vijk( scp.bc().at(b), i,j,k ) { 
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*---------------+
  | edges of walls |
  +---------------*/
  
  /* line i-min & j-min */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=si()-1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = +1;
    xpos = -0.5;
    ypos = -0.5;
    for_k(k) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-min & j-max */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=si()-1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = -1;
    xpos = -0.5;
    ypos = 1.5;
    for_k(k) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-max & j-min */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = +1;
    xpos = 1.5;
    ypos = -0.5;
    for_k(k) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-max & j-max */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = -1;
    xpos = 1.5;
    ypos = 1.5;
    for_k(k) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-min & k-min */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = +1;
    xpos = -0.5;
    zpos = -0.5;
    for_j(j) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-min & k-max */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = -1;
    xpos = -0.5;
    zpos = 1.5;
    for_j(j) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-max & k-min */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofz = +1;
    xpos = 1.5;
    zpos = -0.5;
    for_j(j) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line i-max & k-max */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_j(j) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line j-min & k-min */
  if(scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=sj()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = +1;
    ypos = -0.5;
    zpos = -0.5;
    for_i(i) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line j-min & k-max */
  if(scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=sj()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = -1;
    ypos = -0.5;
    zpos = 1.5;
    for_i(i) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line j-max & k-min */
  if(scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = +1;
    ypos = 1.5;
    zpos = -0.5;
    for_i(i) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /* line j-max & k-max */
  if(scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_i(i) {
      scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }

  /*-----------------+
  | corners of walls |
  +-----------------*/

  /* corner i-min & j-min & k-min */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }


  /* corner i-min & j-min & k-max */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }

  /* corner i-min & j-max & k-min */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }


  /* corner i-min & j-max & k-max */
  if(scp.bc().type_here(Dir::imin(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }


  /* corner i-max & j-min & k-min */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }

  /* corner i-max & j-min & k-max */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmin(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }

  /* corner i-max & j-max & k-min */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }

  /* corner i-max & j-max & k-max */
  if(scp.bc().type_here(Dir::imax(), BndType::wall()) &&
     scp.bc().type_here(Dir::jmax(), BndType::wall()) &&
     scp.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
  }

  /*--------------+
  | immersed body |
  +--------------*/
#if 1
  if(extrapolate_ib) {
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) {

        /* west is in fluid domain */
        if(dom->ibody().on(i-1,j,k)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofx = -1;
          xpos = 1.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

        /* east is in fluid domain */
        if(dom->ibody().on(i+1,j,k)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofx = +1;
          xpos = -0.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

        /* south is in fluid domain */
        if(dom->ibody().on(i,j-1,k)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofy = -1;
          ypos = 1.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

        /* north is in fluid domain */
        if(dom->ibody().on(i,j+1,k)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofy = +1;
          ypos = -0.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

        /* bottom is in fluid domain */
        if(dom->ibody().on(i,j,k-1)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofz = -1;
          zpos = 1.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

        /* top is in fluid domain */
        if(dom->ibody().on(i,j,k+1)) {
          int ofx(0), ofy(0), ofz(0);
          real xpos(0.5), ypos(0.5), zpos(0.5);
          ofz = +1;
          zpos = -0.5;
          scp[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
        }

      } /* center off */
    } /* ijk */
  } /* ibody exists */
  }
#endif
#if 0 /* this is an old version, which doesnt work if
         processor boundary = solid-fluid boundary */
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in solid domain */
    if(dom->ibody().off(i-1,j,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofx = +1;
      xpos = -0.5;
      scp[i-1][j][k] = extrapolate_v(i-1,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofx = -1;
      xpos = 1.5;
      scp[i+1][j][k] = extrapolate_v(i+1,j,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = +1;
      ypos = -0.5;
      scp[i][j-1][k] = extrapolate_v(i,j-1,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = -1;
      ypos = 1.5;
      scp[i][j+1][k] = extrapolate_v(i,j+1,k,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofz = +1;
      zpos = -0.5;
      scp[i][j][k-1] = extrapolate_v(i,j,k-1,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = -1;
      ypos = 1.5;
      scp[i][j][k+1] = extrapolate_v(i,j,k+1,ofx,ofy,ofz,xpos,ypos,zpos,scp);
    }
  }
#endif

  scp.bnd_update_nowall();
  scp.exchange_all();
  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();
 
  return;
}
