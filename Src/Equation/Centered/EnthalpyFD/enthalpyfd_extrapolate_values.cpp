#include "enthalpyfd.h"

/* order of extrapolation */
const int max_ord(3);

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
    if(fin==ctp.idx&&ctp.idx<5)
      cut_stencil.push_back(ctp);

    /* for the purposes of taylor extrapolation, stencil edge is set to zero */
    const real zeropos = cut_stencil[0].pos;
    for(auto & s : cut_stencil)
      s.pos -= zeropos;

    /* point extrapolation */
    for(int idx(0); idx<ctm.idx; ++idx) {
  
      /* extend stencil linearly */
      stencil[idx].pos = stencil[ctm.idx].pos
             - real(ctm.idx-idx)*(stencil[ctm.idx+1].pos-stencil[ctm.idx].pos);
    }
    /* value extrapolation */
    //point_extrapolation(stencil,0,ctm.idx,cut_stencil,zeropos);
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
    if(fin==ctm.idx&&ctm.idx>0)
      cut_stencil.push_back(ctm);

    /* for the purposes of taylor extrapolation, stencil edge is set to zero */
    const real zeropos = cut_stencil[0].pos;
    for(auto & s : cut_stencil)
      s.pos -= zeropos;

    /* point extrapolation */
    for(int idx(ctp.idx+1); idx<=5; ++idx) {

      /* extend stencil linearly */
      stencil[idx].pos = stencil[ctp.idx].pos
             + real(idx-ctp.idx)*(stencil[ctp.idx].pos-stencil[ctp.idx-1].pos);
      /* value extrapolation */
      //stencil[idx].val = point_extrapolation(cut_stencil,
      //                                       stencil[idx].pos-zeropos);
    }
    //point_extrapolation(stencil,ctp.idx+1,5+1,cut_stencil,zeropos);
  }

  return;
}

#if 1
/***************************************************************************//**
*  Lagrangian extrapolation
*******************************************************************************/
real EnthalpyFD::point_extrapolation(const std::vector<StencilPoint> & stencil,
                                     const real xpos) {
  real L(0.0);
  for(int j(0); j<stencil.size(); ++j) {
    real p(1.);
    for(int m(0); m<stencil.size(); ++m) {
      if(m!=j)
        p *= (xpos-stencil[m].pos)/(stencil[j].pos-stencil[m].pos);
    }
    L += stencil[j].val*p;
  }

  return L; 
}
#else
/***************************************************************************//**
*  Taylor extrapolation
*******************************************************************************/
real EnthalpyFD::point_extrapolation(std::vector<StencilPoint> & stencil,
                                     const int i0, const i1,
                                     const std::vector<StencilPoint> & extcil,
                                     const real zeropos) {
  /* construct Taylor polynomial */
  std::vector<real> coefs;
  cht.topo->construct_taylor_polynomial(extcil,coefs);

  /* for linear and quadratic, just extrapolate */
  if(coefs.size()<max_ord+1) {
    for(int idx(i0); idx<i1; ++idx) {
      stencil[idx].val = cht.topo->evaluate_polynomial(coefs.size()-1,coefs,
                                                   stencil[idx].pos-zeropos);
    }
  } else {
    /* test for oscillations */
    real det = 4.*coefs[2]*coefs[2]-12.*coefs[3]*coefs[1];

    /* no extrema or degenerate cubic: just extrapolate */
    if(det<0.||fabs(coefs[3])<boil::atto) {
      for(int idx(i0); idx<i1; ++idx) {
        stencil[idx].val = cht.topo->evaluate_polynomial(coefs.size()-1,coefs,
                                                    stencil[idx].pos-zeropos);
      }
    } else {
      /* write down extrema */
      real ext_1 = (-2.*coefs[2]+sqrt(det))/(6.*coefs[3]);
      real ext_2 = (-2.*coefs[2]-sqrt(det))/(6.*coefs[3]);
      Range<real> ext_range(
                        std::min(extcil->begin().pos,extcil->rbegin().pos),
                        std::max(extcil->begin().pos,extcil->rbegin().pos));

      /* is extcil in an order? */
      bool ext_ordered =  std::is_sorted(extcil.begin(),extcil.end(),
                                         [](const StencilPoint & s1,
                                            const StencilPoint & s2) {
                                              return s1.val<s2.val;  })
                        ||std::is_sorted(extcil.begin(),extcil.end(),
                                         [](const StencilPoint & s1,
                                            const StencilPoint & s2) {
                                              return s1.val>s2.val;  });

      /* if ordered and extremum exists in bounds, go with quadratic */
      if(ext_ordered&&(ext_range.contains(ext_1)||ext_range.contains(ext_2))) {
        for(int idx(i0); idx<i1; ++idx) {
          stencil[idx].val = cht.topo->evaluate_polynomial(coefs.size()-1,coefs,
                                                      stencil[idx].pos-zeropos);
        }
      } else {

    }

  }

#endif
