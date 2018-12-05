#include "finescalar.h"

/******************************************************************************/
real FineScalar::eval_marker(
         const real phi1,
         const int i1, const int j1, const int k1, const int dir1,
         const real phi2,
         const int i2, const int j2, const int k2, const int dir2) {
/***************************************************************************//**
*  \brief Get the marker at the given point.
*******************************************************************************/

  /* degenerate cases */
  if(  (phi1<boil::pico&&phi2-1.0>-boil::pico)
     ||(phi1-1.0>-boil::pico&&phi2<boil::pico)) {
    return 0.5;
  } else {
    point_coord point_1 = eval_point(i1,j1,k1, dir1);
    point_coord point_2 = eval_point(i2,j2,k2, dir2);

    if(point_1.marker == point_2.marker) {
      /* 1 when marker is 1, 0.5 when marker is 0, 0 otherwise */
      return (point_1.marker>0) + 0.5*(point_1.marker==0);
    } else {
      /* they are different */
      if     ( point_1.rl && !point_2.rl) /* 1 is real */
        return point_1.marker>0;
      else if(!point_1.rl &&  point_2.rl) /* 2 is real */
        return point_2.marker>0;
      else if(!point_1.rl && !point_2.rl) /* neither are real */
        return 0.5;
      else
        /* further from interface assumed more valid */
        return (point_1.dist>point_1.dist) ? 
               (point_1.marker>0) :
               (point_2.marker>0) ;
    }
  }
}

/******************************************************************************/
FineScalar::point_coord FineScalar::eval_point(
                 const int i, const int j, const int k, const int dir) {
/***************************************************************************//**
*  \brief Calculate characteristics of the given point.
*******************************************************************************/

  point_coord p;

  /* calculate alpha */
  real alpha = (*nalpha)[i][j][k];
  
  /* degenerate case */
  if(!boil::realistic(alpha)) {
    real c = (*phi)[i][j][k];
    if(c>phisurf)
     p.marker = 1;
    else
     p.marker = -1;

    p.dist = 0.0;
    p.rl = false;

    return p;
  }

  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -(*nx)[i][j][k];
  real vn2 = -(*ny)[i][j][k];
  real vn3 = -(*nz)[i][j][k];

  /* calculate position of interest */
  real xpos(0.5),ypos(0.5),zpos(0.5);

  if(dir == ws()) { xpos = 0.0; ypos = 0.0; }
  if(dir == wn()) { xpos = 0.0; ypos = 1.0; }
  if(dir == es()) { xpos = 1.0; ypos = 0.0; }
  if(dir == en()) { xpos = 1.0; ypos = 1.0; }

  if(dir == wb()) { xpos = 0.0; zpos = 0.0; }
  if(dir == wt()) { xpos = 0.0; zpos = 1.0; }
  if(dir == eb()) { xpos = 1.0; zpos = 0.0; }
  if(dir == et()) { xpos = 1.0; zpos = 1.0; }
        
  if(dir == sb()) { ypos = 0.0; zpos = 0.0; }
  if(dir == st()) { ypos = 0.0; zpos = 1.0; }
  if(dir == nb()) { ypos = 1.0; zpos = 0.0; }
  if(dir == nt()) { ypos = 1.0; zpos = 1.0; }

  if(dir == wsb()) { xpos = 0.0; ypos = 0.0; zpos = 0.0; }
  if(dir == wst()) { xpos = 0.0; ypos = 0.0; zpos = 1.0; }
  if(dir == wnb()) { xpos = 0.0; ypos = 1.0; zpos = 0.0; }
  if(dir == wnt()) { xpos = 0.0; ypos = 1.0; zpos = 1.0; }
  if(dir == esb()) { xpos = 1.0; ypos = 0.0; zpos = 0.0; }
  if(dir == est()) { xpos = 1.0; ypos = 0.0; zpos = 1.0; }
  if(dir == enb()) { xpos = 1.0; ypos = 1.0; zpos = 0.0; }
  if(dir == ent()) { xpos = 1.0; ypos = 1.0; zpos = 1.0; }


  /* mirror position of interest */
  if(vn1<0) xpos = 1.0-xpos;
  if(vn2<0) ypos = 1.0-ypos;
  if(vn3<0) zpos = 1.0-zpos;

  /* calculate mirrored normal vector */
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  /* calculate distance to the interface */
  real dist = alpha - vm1*xpos - vm2*ypos - vm3*zpos;

  /* decide on position wrt interface */
  /* vm points to gas */
  int marker(0);
  if(dist>boil::pico) {
    marker = 1;  /* liquid */
  } else if(dist<-boil::pico) {
    marker = -1; /* gas */
    dist = fabs(dist);
  } else {
    dist = 0.0;
  }

  p.marker = marker;
  p.dist = dist;

  real c = (*phi)[i][j][k];
  /* is there interface bw centre and point? */
  if((c>=phisurf&&marker<0)||(c<=phisurf&&marker>0))
    p.rl = true;
  else
    p.rl = false;

  return p;

}
