#include "enthalpyfd.h"

#include <random>
#include <iomanip>

/******************************************************************************/
bool EnthalpyFD::test_extrapolation(const int count) {
/******************************************************************************/

  std::mt19937 rng(123);
  std::uniform_real_distribution<real> distR1(1.,10.);
  std::uniform_real_distribution<real> distR2(-10.,10.);
  std::uniform_int_distribution<int> distI1(0,1);
  std::uniform_int_distribution<int> distI2(0,3);
  std::uniform_int_distribution<int> distI3(2,5);

  bool flag(true);
  for(int i(0); i<count; ++i) {

    int i1 = distI2(rng);
    int i2 = distI3(rng);

    boil::oout<<"--------------- "<<i<<" | "<<i1<<" "<<i2<<boil::endl;

    std::vector<real> stencil = { (distI1(rng)>0 ? 1 : -1) * distR1(rng),
                                  (distI1(rng)>0 ? 1 : -1) * distR1(rng),
                                  (distI1(rng)>0 ? 1 : -1) * distR1(rng),
                                  (distI1(rng)>0 ? 1 : -1) * distR1(rng),
                                  (distI1(rng)>0 ? 1 : -1) * distR1(rng),
                                  (distI1(rng)>0 ? 1 : -1) * distR1(rng)};

    /* completely randomised stencil with random cubic polynomial */
    flag = flag & test_extrapolation(stencil,
                                     {distR2(rng),distR2(rng),
                                      distR2(rng),distR2(rng)},
                                     {std::min(i1,i2), std::max(i1,i2)});
  }

  return flag;

}

/******************************************************************************/
bool EnthalpyFD::test_extrapolation(std::vector<real> & stencil,
                                    const std::vector<real> & coefficients,
                                    const std::vector<int> & cutpoints) {
/***************************************************************************//**
*  \brief Test if extrapolation works correctly using polynomials
*******************************************************************************/
  bool flag(true);

  /* sort stencil */
  std::sort(stencil.begin(),stencil.end());

  //std::cout.setf(std::ios_base::scientific);
  std::cout<< std::setprecision(3);

  /* output */
  boil::oout<<"S ";
  for(auto s : stencil)
    boil::oout<<s<<" ";
  boil::oout<<"| C ";
  for(auto c : coefficients)
    boil::oout<<c<<" ";
  boil::oout<<"| I ";
  for(auto c : cutpoints)
    boil::oout<<c<<" ";
  boil::oout<<boil::endl;

  /* we construct a polynomial of given order. 
   * Then we construct a Lagrangian polynomial of the same order and
   * extrapolate to another point. If this works, values should be 
   * exactly the same */

  int ord(3);
  if(cutpoints[0]==cutpoints[1])
    ord = 2;
  
  /* fill stencil */
  std::vector<StencilPoint> stp,sttest;
  int i(0);
  for(auto & s : stencil) {
    real val = cht.topo->evaluate_polynomial(ord,coefficients,s);
    stp.push_back(StencilPoint(i,val,s));
    ++i;
  }
  
  /* create cutoffs */
  real ctm_pos = (cutpoints[0]>0) 
               ? 0.5*(stencil[cutpoints[0]-1]+stencil[cutpoints[0]])
               : stencil[0];
  real ctp_pos = (cutpoints[1]<5) 
               ? 0.5*(stencil[cutpoints[1]]+stencil[cutpoints[1]+1])
               : stencil[5];

  real ctm_val = cht.topo->evaluate_polynomial(ord,coefficients,ctm_pos);
  real ctp_val = cht.topo->evaluate_polynomial(ord,coefficients,ctp_pos);

  StencilPoint ctm(cutpoints[0],ctm_val,ctm_pos), 
               ctp(cutpoints[1],ctp_val,ctp_pos);

  boil::oout<<"Cutoffs: "<<ctm.idx<<" "<<ctm.pos<<" "<<ctm.val
            <<" "<<cht.topo->evaluate_polynomial_derivative(ord,coefficients,ctm.pos)
            <<" | ";
  boil::oout<<ctp.idx<<" "<<ctp.pos<<" "<<ctp.val
            <<" "<<cht.topo->evaluate_polynomial_derivative(ord,coefficients,ctp.pos)
            <<" | "<<boil::endl;

  /* test */
  boil::oout<<"pol-ders | ";
  for(auto & s : stp) {
    boil::oout<<s.idx<<" "<<s.pos
              <<" "<<cht.topo->evaluate_polynomial_derivative(ord,coefficients,s.pos)
              <<" | ";
  }
  boil::oout<<boil::endl;

  boil::oout<<"pre-test | ";
  for(auto & s : stp) {
    boil::oout<<s.idx<<" "<<s.pos<<" "<<s.val
              <<" | ";
  }
  boil::oout<<boil::endl;

  extrapolate_values(stp,ctm,ctp);

  /* control */
  i = 0;
  for(auto & s : stp) {
    real val = cht.topo->evaluate_polynomial(ord,coefficients,s.pos);
    sttest.push_back(StencilPoint(i,val,s.pos));
    ++i;
  }

  boil::oout<<"post-test| ";
  for(auto & s : stp) {
    boil::oout<<s.idx<<" "<<s.pos<<" "<<s.val<<" | ";
  }
  boil::oout<<boil::endl;
  boil::oout<<"control  | ";
  for(auto & s : sttest) {
    boil::oout<<s.idx<<" "<<s.pos<<" "<<s.val<<" | ";
  }
  boil::oout<<boil::endl;

  /* evaluate */
  for(int idx(0); idx<sttest.size();++idx) {
    flag &= (fabs(sttest[idx].val-stp[idx].val)
                  /std::max(fabs(sttest[idx].val),boil::nano))<boil::micro;

    if(!flag) {
      boil::oout<<"Failed: "<<fabs(sttest[idx].val-stp[idx].val)
                <<" "<<(fabs(sttest[idx].val-stp[idx].val)
                       /std::max(fabs(sttest[idx].val),boil::nano))
                <<" "<<sttest[idx].val<<" "<<stp[idx].val<<boil::endl;
      exit(0);
    }
  }

  return flag;
}
