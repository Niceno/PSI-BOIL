#include "enthalpyfd.h"

/***************************************************************************//**
* Calculate the value of phi at faces using a high-order non-linear scheme.
*******************************************************************************/
real EnthalpyFD::face_value(const Sign matter_sig, const Comp m, const real vel,
                            const int i, const int j, const int k,
                            const int ofx, const int ofy, const int ofz,
                            const Old old) {

  /* dummy sign */
  Sign dummy;
   
  /* are we in a liquid-cell or near-liquid cell? */
  bool west_above(true), east_above(true);
  if(cht.above_interface(i-ofx,j-ofy,k-ofz,old)) {
    if(!cht.above_interface(i,j,k,old)) {
      east_above = false;
    }
  } else if(cht.above_interface(i,j,k,old)) {
    west_above = false;
  } else {
    west_above = east_above = false;
  }

  /* do we solve for gas? -> inversion */
  if(matter_sig<0) {
    west_above = !west_above;
    east_above = !east_above;
  }

  /* are we in a region, in which we don't need the respective flux? */
  if(!west_above&&!east_above) 
    return 0.0;

  /* edge of domain */
  if(stencil_min(m,i,j,k)) {
    return phi[i-ofx][j-ofy][k-ofz];
  }
  if(stencil_max(m,i-ofx,j-ofy,k-ofz)) {
    return phi[i][j][k];
  }

  /* ibody */
  if(dom->ibody().off(i,j,k)||dom->ibody().off(i-ofx,j-ofy,k-ofz)) {
    return 0.0;
  }

  /* 6-stencil, offset due to staggering */
  std::vector<StencilPoint> stencil;
  for(int s(0); s<6; ++s)
    stencil.push_back(StencilPoint(s));

  stencil[2].pos = -cht.distance_face(Sign::pos(),m,i-ofx,j-ofy,k-ofz);
  stencil[1].pos = stencil[2].pos
            - cht.distance_center(Sign::neg(),m,i-ofx,j-ofy,k-ofz);
  stencil[0].pos = stencil[1].pos
            - cht.distance_center(Sign::neg(),m,i-2*ofx,j-2*ofy,k-2*ofz);
  stencil[3].pos = cht.distance_face(Sign::neg(),m,i,j,k);
  stencil[4].pos = stencil[3].pos
            + cht.distance_center(Sign::pos(),m,i,j,k);
  stencil[5].pos = stencil[4].pos
            + cht.distance_center(Sign::pos(),m,i+ofx,j+ofy,k+ofz);

  stencil[0].val = phi[i-3*ofx][j-3*ofy][k-3*ofz];
  stencil[1].val = phi[i-2*ofx][j-2*ofy][k-2*ofz];
  stencil[2].val = phi[i  -ofx][j  -ofy][k  -ofz];
  stencil[3].val = phi[i      ][j      ][k      ];
  stencil[4].val = phi[i  +ofx][j  +ofy][k  +ofz];
  stencil[5].val = phi[i+2*ofx][j+2*ofy][k+2*ofz];

  /* stencil cutoff: lowest/highest correct idx, followed by the cutoff */
  StencilPoint ctm(0),ctp(5);

  /* west */
  if(west_above) {
    /* m3 */
    if(cht.interface(Sign::neg(),m,i-2*ofx,j-2*ofy,k-2*ofz,old)) {
      ctm = StencilPoint(1);
      ctm.pos = stencil[1].pos
              - cht.distance_int(Sign::neg(),m,i-2*ofx,j-2*ofy,k-2*ofz,
                                 ctm.val,dummy,ResistEval::yes,old);
    } else if(stencil_min(m,i-2*ofx,j-2*ofy,k-2*ofz)) {
      ctm = StencilPoint(1);
      ctm.pos = stencil[0].pos;
      ctm.val = stencil[0].val;
    } else if(dom->ibody().off(i-3*ofx,j-3*ofy,k-3*ofz)) {
      ctm = StencilPoint(1);
      ctm.pos = stencil[1].pos
              - cht.distance_face(Sign::neg(),m,i-2*ofx,j-2*ofy,k-2*ofz);
      ctm.val = cht.node_tmp_flu()[m][i-2*ofx][j-2*ofy][k-2*ofz]; 
    }

    /* m2 */
    if(cht.interface(Sign::neg(),m,i-ofx,j-ofy,k-ofz,old)) {
      ctm = StencilPoint(2);
      ctm.pos = stencil[2].pos
              - cht.distance_int(Sign::neg(),m,i-ofx,j-ofy,k-ofz,
                                 ctm.val,dummy,ResistEval::yes,old);
    } else if(stencil_min(m,i-ofx,j-ofy,k-ofz)) {
      ctm = StencilPoint(2);
      ctm.pos = stencil[1].pos;
      ctm.val = stencil[1].val;
    } else if(dom->ibody().off(i-2*ofx,j-2*ofy,k-2*ofz)) {
      ctm = StencilPoint(2);
      ctm.pos = stencil[2].pos
             - cht.distance_face(Sign::neg(),m,i-ofx,j-ofy,k-ofz);
      ctm.val = cht.node_tmp_flu()[m][i-ofx][j-ofy][k-ofz];
    }

  } else { /* west above? */

    /* m1 */
    if(cht.interface(Sign::neg(),m,i,j,k,old)) {
      ctm = StencilPoint(3);
      ctm.pos = stencil[3].pos
              - cht.distance_int(Sign::neg(),m,i,j,k,
                                 ctm.val,dummy,ResistEval::yes,old);
    }

  }

  /* east */
  if(east_above) {
    /* p3 */
    if(cht.interface(Sign::pos(),m,i+ofx,j+ofy,k+ofz,old)) {
      ctp = StencilPoint(4);
      ctp.pos = stencil[4].pos
              + cht.distance_int(Sign::pos(),m,i+ofx,j+ofy,k+ofz,
                                 ctp.val,dummy,ResistEval::yes,old);
    } else if(stencil_max(m,i+ofx,j+ofy,k+ofz)) {
      ctp = StencilPoint(4);
      ctp.pos = stencil[5].pos;
      ctp.val = stencil[5].val;
    } else if(dom->ibody().off(i+2*ofx,j+2*ofy,k+2*ofz)) {
      ctp = StencilPoint(4);
      ctp.pos = stencil[4].pos
             + cht.distance_face(Sign::pos(),m,i+ofx,j+ofy,k+ofz);
      ctp.val = cht.node_tmp_flu()[m][i+2*ofx][j+2*ofy][k+2*ofz];
    }

    /* p2 */
    if(cht.interface(Sign::pos(),m,i,j,k,old)) {
      ctp = StencilPoint(3);
      ctp.pos = stencil[3].pos
              + cht.distance_int(Sign::pos(),m,i,j,k,
                                 ctp.val,dummy,ResistEval::yes,old);
    } else if(stencil_max(m,i,j,k)) {
      ctp = StencilPoint(3);
      ctp.pos = stencil[4].pos;
      ctp.val = stencil[4].val;
    } else if(dom->ibody().off(i+ofx,j+ofy,k+ofz)) {
      ctp = StencilPoint(3);
      ctp.pos = stencil[3].pos
              + cht.distance_face(Sign::pos(),m,i,j,k);
      ctp.val = cht.node_tmp_flu()[m][i+ofx][j+ofy][k+ofz];
    }

  } else { /* east above? */

    /* p1 */
    if(cht.interface(Sign::pos(),m,i-ofx,j-ofy,k-ofz,old)) {
      ctp = StencilPoint(2);
      ctp.pos = stencil[2].pos
              + cht.distance_int(Sign::pos(),m,i-ofx,j-ofy,k-ofz,
                                 ctp.val,dummy,ResistEval::yes,old);
    }
  }

  /* sanity check */
  assert(ctp.idx>=ctm.idx);

#if 1
  /* extrapolation of missing values */
  extrapolate_values(stencil,ctm,ctp);

  /* calculate face value of temperature */
  real tval;

  /* upwind TVD */
  if(vel>0.) {
    tval = lim.limit(+1., stencil[1].val, stencil[2].val, stencil[3].val);
  } else {
    tval = lim.limit(+1., stencil[4].val, stencil[3].val, stencil[2].val);
  }

  /* evaluate flux */
  return tval;
#else
  /* upwind TVD */
  if(vel>0.) {
    if(ctm.idx<=1&&ctp.idx>=3) {
  #if 1
      return lim.limit(+1., stencil[1].val, stencil[2].val, stencil[3].val);
  #else
      ctp = StencilPoint(2,stencil[3].val,stencil[3].pos);
      return extrapolate_value(Sign::pos(),stencil,ctm,ctp,0.);
  #endif
    } else if(fabs(real(ctm.idx)-2.5)>fabs(real(ctp.idx)-2.5)) {
      return extrapolate_value(Sign::pos(),stencil,ctm,ctp,0.);
    } else {
      return extrapolate_value(Sign::neg(),stencil,ctm,ctp,0.);
    }
  } else {
    if(ctp.idx>=4&&ctm.idx<=2) {
  #if 1
      return lim.limit(+1., stencil[4].val, stencil[3].val, stencil[2].val);
  #else
      ctm = StencilPoint(3,stencil[2].val,stencil[2].pos);
      return extrapolate_value(Sign::neg(),stencil,ctm,ctp,0.);
  #endif
    } else if(fabs(real(ctp.idx)-2.5)>fabs(real(ctm.idx)-2.5)) {
      return extrapolate_value(Sign::neg(),stencil,ctm,ctp,0.);
    } else {
      return extrapolate_value(Sign::pos(),stencil,ctm,ctp,0.);
    }
  }

  return 0.0;
#endif
}
