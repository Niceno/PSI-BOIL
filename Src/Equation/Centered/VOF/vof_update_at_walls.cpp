#include "vof.h"

/******************************************************************************/
void VOF::update_at_walls() {
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
  for( int b=0; b<phi.bc().count(); b++ ) {

    if( phi.bc().type_decomp(b) ) continue;

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
          phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
        }
      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */

  /*---------------+
  | edges of walls |
  +---------------*/
  
  /* line i-min & j-min */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=si()-1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = +1;
    xpos = -0.5;
    ypos = -0.5;
    for_k(k) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & j-max */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=si()-1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofy = -1;
    xpos = -0.5;
    ypos = 1.5;
    for_k(k) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & j-min */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei()+1;
    int j=sj()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = +1;
    xpos = 1.5;
    ypos = -0.5;
    for_k(k) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & j-max */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei()+1;
    int j=ej()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofy = -1;
    xpos = 1.5;
    ypos = 1.5;
    for_k(k) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & k-min */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = +1;
    xpos = -0.5;
    zpos = -0.5;
    for_j(j) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-min & k-max */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = +1;
    ofz = -1;
    xpos = -0.5;
    zpos = 1.5;
    for_j(j) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & k-min */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofx = -1;
    ofz = +1;
    xpos = 1.5;
    zpos = -0.5;
    for_j(j) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line i-max & k-max */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_j(j) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-min & k-min */
  if(phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=sj()-1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = +1;
    ypos = -0.5;
    zpos = -0.5;
    for_i(i) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-min & k-max */
  if(phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=sj()-1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = +1;
    ofz = -1;
    ypos = -0.5;
    zpos = 1.5;
    for_i(i) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-max & k-min */
  if(phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej()+1;
    int k=sk()-1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = +1;
    ypos = 1.5;
    zpos = -0.5;
    for_i(i) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /* line j-max & k-max */
  if(phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej()+1;
    int k=ek()+1;
    int ofx(0), ofy(0), ofz(0);
    real xpos(0.5), ypos(0.5), zpos(0.5);
    ofy = -1;
    ofz = -1;
    ypos = 1.5;
    zpos = 1.5;
    for_i(i) {
      phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  /*-----------------+
  | corners of walls |
  +-----------------*/

  /* corner i-min & j-min & k-min */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-min & j-min & k-max */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-min & j-max & k-min */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-min & j-max & k-max */
  if(phi.bc().type_here(Dir::imin(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }


  /* corner i-max & j-min & k-min */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-min & k-max */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmin(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-max & k-min */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmin(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

  /* corner i-max & j-max & k-max */
  if(phi.bc().type_here(Dir::imax(), BndType::wall()) &&
     phi.bc().type_here(Dir::jmax(), BndType::wall()) &&
     phi.bc().type_here(Dir::kmax(), BndType::wall())   ) {
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
    phi[i][j][k] = extrapolate_v(i,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
  }

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
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofx = +1;
      xpos = -0.5;
      phi[i-1][j][k] = extrapolate_v(i-1,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }

    /* east */
    if (dom->ibody().off(i+1,j,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofx = -1;
      xpos = 1.5;
      phi[i+1][j][k] = extrapolate_v(i+1,j,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }

    /* south */
    if (dom->ibody().off(i,j-1,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = +1;
      ypos = -0.5;
      phi[i][j-1][k] = extrapolate_v(i,j-1,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }

    /* north */
    if (dom->ibody().off(i,j+1,k)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = -1;
      ypos = 1.5;
      phi[i][j+1][k] = extrapolate_v(i,j+1,k,ofx,ofy,ofz,xpos,ypos,zpos);
    }

    /* bottom */
    if (dom->ibody().off(i,j,k-1)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofz = +1;
      zpos = -0.5;
      phi[i][j][k-1] = extrapolate_v(i,j,k-1,ofx,ofy,ofz,xpos,ypos,zpos);
    }

    /* top */
    if (dom->ibody().off(i,j,k+1)) {
      int ofx(0), ofy(0), ofz(0);
      real xpos(0.5), ypos(0.5), zpos(0.5);
      ofy = -1;
      ypos = 1.5;
      phi[i][j][k+1] = extrapolate_v(i,j,k+1,ofx,ofy,ofz,xpos,ypos,zpos);
    }
  }

  phi.exchange_all();
  //nx.exchange_all();
  //ny.exchange_all();
  //nz.exchange_all();
 
  return;
}

/*-------------------+
| ancillary function |
+-------------------*/
real VOF::extrapolate_v(const int i, const int j, const int k,
                        const int ofx, const int ofy, const int ofz,
                        const real xp, const real yp, const real zp) {

  int ii = i+ofx;
  int jj = j+ofy;
  int kk = k+ofz;
   
  real phiphi = phi[ii][jj][kk];

  /* erroneous interfaces */
  if(phiphi<tol_wall||phiphi-1.0>-tol_wall) {
    return real(phiphi>phisurf);
  } 

  /* unnormalized alpha value */
  real alphaval = nalpha[ii][jj][kk];
          
  /* degenerate case I */
  if(!boil::realistic(alphaval)) {
    return real(phiphi>phisurf);
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
    return real(phiphi>phisurf);
  }
         
  real xpos = xp;
  real ypos = yp;
  real zpos = zp;

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
  vm1 /= denom;
  vm2 /= denom;
  vm3 /= denom;

#if 0
  if(i==1&&j==0&&k==90) boil::oout<<ii<<" "<<jj<<" "<<kk<<" "<<alphaval<<" "<<vm1*(0.5-xpos) + vm2*(0.5-ypos) + vm3*(0.5-zpos)<<" "<<phiphi<<" "<<calc_v(alphaval,vm1,vm2,vm3)<<" | "<<vn1<<" "<<vn2<<" "<<vn3<<" "<<xpos<<" "<<ypos<<" "<<zpos<<" | "<<ofx<<" "<<ofy<<" "<<ofz<<boil::endl;
  //if(i==0&&j==1&&k==89) boil::oout<<ii<<" "<<jj<<" "<<kk<<" "<<alphaval<<" "<<vm1*(0.5-xpos) + vm2*(0.5-ypos) + vm3*(0.5-zpos)<<" "<<phiphi<<" "<<calc_v(alphaval,vm1,vm2,vm3)<<" | "<<vn1<<" "<<vn2<<" "<<vn3<<" "<<xpos<<" "<<ypos<<" "<<zpos<<" | "<<ofx<<" "<<ofy<<" "<<ofz<<boil::endl;
#endif

  /* volume fraction */
  return calc_v(alphaval,vm1,vm2,vm3);

}
