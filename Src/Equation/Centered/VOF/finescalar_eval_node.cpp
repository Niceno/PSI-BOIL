#include "finescalar.h"

/******************************************************************************/
void FineScalar::eval_node() {
/***************************************************************************//**
*  \brief Evaluate the marker function at nodes.
*******************************************************************************/

  /* staggered in x and y and z */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mmm = (*phi)[i-1][j-1][k-1];
    real phi_mmp = (*phi)[i-1][j-1][k  ];
    real phi_mpm = (*phi)[i-1][j  ][k-1];
    real phi_mpp = (*phi)[i-1][j  ][k  ];
    real phi_pmm = (*phi)[i  ][j-1][k-1];
    real phi_pmp = (*phi)[i  ][j-1][k  ];
    real phi_ppm = (*phi)[i  ][j  ][k-1];
    real phi_ppp = (*phi)[i  ][j  ][k  ];

    /* evaluate each cell pair */
    bool mmm_ppp = (phi_mmm-phisurf)*(phi_ppp-phisurf)>0.0;
    bool mmp_ppm = (phi_mmp-phisurf)*(phi_ppm-phisurf)>0.0;
    bool mpp_pmm = (phi_mpp-phisurf)*(phi_pmm-phisurf)>0.0;
    bool mpm_pmp = (phi_mpm-phisurf)*(phi_pmp-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[4];
    int idx(0);

    if(mmm_ppp) {
      vals[idx] = phi_mmm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mmm, i-1,j-1,k-1, ent(),
                              phi_ppp, i  ,j  ,k  , wsb());
    }
    idx++;
    if(mmp_ppm) {
      vals[idx] = phi_mmp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mmp, i-1,j-1,k  , enb(),
                              phi_ppm, i  ,j  ,k-1, wst());
    }
    idx++;
    if(mpp_pmm) {
      vals[idx] = phi_mpp>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mpp, i-1,j  ,k  , esb(),
                              phi_pmm, i  ,j-1,k-1, wnt());
    }
    idx++;
    if(mpm_pmp) {
      vals[idx] = phi_mpm>phisurf;
    } else {
      vals[idx] = eval_marker(phi_mpm, i-1,j  ,k-1, est(),
                              phi_pmp, i  ,j-1,k  , wnb());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k, wsb()) = sum/4.0;
  }

  /* correct at boundaries */
  //bdcond_node();

  return;
}

#if 0
/***********************
 * ancillary functions
 ***********************/
real FineScalar::edge_eval_val(
         const real phi1,
         const int i1, const int j1, const int k1, const int dir1,
         const real phi2,
         const int i2, const int j2, const int k2, const int dir2) {
  /* degenerate cases */
  if(  (phi1<boil::pico&&phi2-1.0>-boil::pico)
     ||(phi1-1.0>-boil::pico&&phi2<boil::pico)) {
    return 0.5;
  } else {
    point_coord point_1 = edge_eval_point(i1,j1,k1, dir1);
    point_coord point_2 = edge_eval_point(i2,j2,k2, dir2);

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

FineScalar::point_coord FineScalar::edge_eval_point(
                 const int i, const int j, const int k, const int dir) {

  point_coord p;

  /* calculate alpha */
  real alpha = (*nalpha)[i][j][k];
  
  /* degenerate case */
  if(alpha>boil::zetta) {
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
#endif
