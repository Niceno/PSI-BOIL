#include "finescalar.h"

/******************************************************************************/
void FineScalar::eval_edge() {
/***************************************************************************//**
*  \brief Evaluate the marker function at edges.
*******************************************************************************/

  /* staggered in y and z */
  for(int i=phi->si(); i<=phi->ei()  ; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i][j-1][k-1];
    real phi_mp = (*phi)[i][j-1][k  ];
    real phi_pm = (*phi)[i][j  ][k-1];
    real phi_pp = (*phi)[i][j  ][k  ];

    /* evaluate each cell pair */
#if 0
    bool mm_mp = (phi_mm-phisurf)*(phi_mp-phisurf)>0.0;
    bool mm_pm = (phi_mm-phisurf)*(phi_pm-phisurf)>0.0;
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;
    bool mp_pp = (phi_mp-phisurf)*(phi_pp-phisurf)>0.0;
    bool pm_pp = (phi_pm-phisurf)*(phi_pp-phisurf)>0.0;
#else
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;
#endif

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mm, i,j-1,k-1, nt(),
                                phi_pp, i,j  ,k  , sb());
    }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mp, i,j-1,k  , nb(),
                                phi_pm, i,j  ,k-1, st());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,sb()) = sum/2.0;
     
  }

  /* staggered in x and z */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()  ; j++)
  for(int k=phi->sk(); k<=phi->ek()+1; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i-1][j][k-1];
    real phi_mp = (*phi)[i-1][j][k  ];
    real phi_pm = (*phi)[i  ][j][k-1];
    real phi_pp = (*phi)[i  ][j][k  ];

    /* evaluate each cell pair */
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mm, i-1,j,k-1, et(),
                                phi_pp, i  ,j,k  , wb());
    }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mp, i-1,j,k  , eb(),
                                phi_pm, i  ,j,k-1, wt());
    }

    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,wb()) = sum/2.0;
     
  }

  /* staggered in x and y */
  for(int i=phi->si(); i<=phi->ei()+1; i++)
  for(int j=phi->sj(); j<=phi->ej()+1; j++)
  for(int k=phi->sk(); k<=phi->ek()  ; k++) {

    /* declare/initialize necessary values */
    real phi_mm = (*phi)[i-1][j-1][k];
    real phi_mp = (*phi)[i-1][j  ][k];
    real phi_pm = (*phi)[i  ][j-1][k];
    real phi_pp = (*phi)[i  ][j  ][k];

    /* evaluate each cell pair */
    bool mm_pp = (phi_mm-phisurf)*(phi_pp-phisurf)>0.0;
    bool mp_pm = (phi_mp-phisurf)*(phi_pm-phisurf)>0.0;

    /* evaluate point candidates in question */
    real vals[2];
    int idx(0);

    if(mm_pp) {
      vals[idx] = phi_mm>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mm, i-1,j-1,k, en(),
                                phi_pp, i  ,j  ,k, ws());
     }
    idx++;
    if(mp_pm) {
      vals[idx] = phi_mp>phisurf;
    } else {
      vals[idx] = edge_eval_val(phi_mp, i-1,j  ,k, es(),
                                phi_pm, i  ,j-1,k, wn());
    }
 
    /* sum everything together */
    real sum(0.0);
    for(auto val : vals)
      sum += val;
    value(i,j,k,ws()) = sum/2.0;
     
  }

  /* correct at boundaries */
  //bdcond_edge();

  return;
}

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
