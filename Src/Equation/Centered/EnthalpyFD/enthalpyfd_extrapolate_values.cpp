#include "enthalpyfd.h"

/* order of extrapolation */
const int max_ord(1);

/***************************************************************************//**
*  \brief Extrapolate value stencil for gradient calculation.
*******************************************************************************/
void EnthalpyFD::extrapolate_values(std::vector<StencilPoint> & stencil,
                                    const StencilPoint & ctm,
                                    const StencilPoint & ctp) {

  /* extrapolate west */
  if(ctm.idx>0) {
  
    /* construct the extrapolation stencil */
    std::vector<StencilPoint> cut_stencil;
    int fin = std::min(ctm.idx+max_ord-1,ctp.idx);

    cut_stencil.push_back(ctm);
    for(int idx(ctm.idx); idx<=fin; ++idx) {
      cut_stencil.push_back(stencil[idx]);
    }
    if(cut_stencil.size()<max_ord+1&&ctp.idx<5)
      cut_stencil.push_back(ctp);

    /* point extrapolation */
    for(int idx(0); idx<ctm.idx; ++idx) {
  
      /* extend stencil linearly */
      stencil[idx].pos = stencil[ctm.idx].pos
             - real(ctm.idx-idx)*(stencil[ctm.idx+1].pos-stencil[ctm.idx].pos);
    }
    /* value extrapolation */
    point_extrapolation(stencil,0,ctm.idx,cut_stencil);
  }

  /* extrapolate east */
  if(ctp.idx<5) {

    /* construct the extrapolation stencil */
    std::vector<StencilPoint> cut_stencil;
    int fin = std::max(ctp.idx-max_ord+1,ctm.idx);

    cut_stencil.push_back(ctp);
    for(int idx(ctp.idx); idx>=fin; --idx) {
      cut_stencil.push_back(stencil[idx]);
    }
    if(cut_stencil.size()<max_ord+1&&ctm.idx>0)
      cut_stencil.push_back(ctm);

    /* point extrapolation */
    for(int idx(ctp.idx+1); idx<6; ++idx) {

      /* extend stencil linearly */
      stencil[idx].pos = stencil[ctp.idx].pos
             + real(idx-ctp.idx)*(stencil[ctp.idx].pos-stencil[ctp.idx-1].pos);
    }
    point_extrapolation(stencil,ctp.idx+1,5+1,cut_stencil);
  }

  return;
}

#if 1
/***************************************************************************//**
*  Lagrangian extrapolation
*******************************************************************************/
void EnthalpyFD::point_extrapolation(std::vector<StencilPoint> & stencil,
                                     const int i0, const int i1,
                                     const std::vector<StencilPoint> & extcil) {

  assert(extcil.size()==max_ord+1);

  for(int idx(i0); idx<i1; ++idx) {
    real L(0.0);
    for(int j(0); j<extcil.size(); ++j) {
      real p(1.);
      for(int m(0); m<extcil.size(); ++m) {
        if(m!=j)
          p *= (stencil[idx].pos-extcil[m].pos)/(extcil[j].pos-extcil[m].pos);
      }
      L += extcil[j].val*p;
    }
    stencil[idx].val = L;
  }
}
#else
/***************************************************************************//**
*  Taylor extrapolation. We use the fact that for a polynomial constructed with
*  as many points as the order is, the difference at a given point is exactly
*  equal to the derivative at that point.
*******************************************************************************/
void EnthalpyFD::point_extrapolation(std::vector<StencilPoint> & stencil,
                                     const int i0, const int i1,
                                     const std::vector<StencilPoint> & extcil) {

  std::vector<StencilPoint> modcil = extcil;
  int ord(extcil.size()-1);

  /*** testing for reduced-order ***/
  bool flag =   (ord<=2) /* for linear and quadratic, just extrapolate */
              || /* non-monotonic data */
               !(  std::is_sorted(modcil.begin(),modcil.end(),
                                  [](const StencilPoint & s1,
                                     const StencilPoint & s2) {
                                     return s1.val<s2.val;  })
                 ||std::is_sorted(modcil.begin(),modcil.end(),
                                  [](const StencilPoint & s1,
                                     const StencilPoint & s2) {
                                     return s1.val>s2.val;  })
                );

  if(flag) {
    ord = std::min(2,ord);
  } else {
    /* for monotonic, we test the derivatives at stencil point */
    std::vector<real> ders;

    for(int idx(0); idx<modcil.size(); ++idx) {
      /* shift origin to extrapolated point */
      for(int jdx(0); jdx<modcil.size(); ++jdx)
        modcil[jdx].pos = extcil[jdx].pos-extcil[idx].pos;

      /* evaluate */
      ders.push_back(cht.topo->nth_order_first(modcil,
                                               AccuracyOrder(ord)));
    }
    
    /* is the derivative of the same sign? */
    bool mono(true);
    for(int idx(0); idx<ders.size()-1; ++idx) {
      mono &= signum(1.0,ders[idx])==signum(1.0,ders[idx+1]);
    }

    /* fall-back to quadratic */
    if(!mono)
      ord = std::min(2,ord);

  } /* if flag */

  /* extrapolate */
  for(int idx(i0); idx<i1; ++idx) {
    /* shift origin to extrapolated point */
    for(int jdx(0); jdx<modcil.size(); ++jdx) {
      modcil[jdx].pos = extcil[jdx].pos-stencil[idx].pos;
    }

    /* extrapolate */
    stencil[idx].val = cht.topo->nth_order_zeroth(modcil,
                                                  AccuracyOrder(ord));
  }

  return;
}
#endif
