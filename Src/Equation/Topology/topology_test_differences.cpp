#include "topology.h"

#include <random>

/******************************************************************************/
bool Topology::test_differences_first(const int count) {
/******************************************************************************/
  std::mt19937 rng(123);
  std::uniform_real_distribution<real> distR1(1.,10.);
  std::uniform_real_distribution<real> distR2(-10.,10.);
  std::uniform_int_distribution<int> distI(0,1);

  bool flag(true);
  for(int i(0); i<count; ++i) {
 
    boil::oout<<i<<boil::endl;

    /* completely randomised stencil with random polynomial */
    flag = flag & test_differences_first(
                                   {0.,(distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng)},
                                   {distR2(rng),distR2(rng),distR2(rng),
                                    distR2(rng),distR2(rng)});
  }

  return flag;

}


/******************************************************************************/
bool Topology::test_differences_first(const std::vector<real> & stenpos,
                                      const std::vector<real> & coefficients) {
/***************************************************************************//**
*  \brief Test if differences work correctly using polynomials
*******************************************************************************/

  bool flag(true);

  for(auto s : stenpos)
    boil::oout<<s<<" ";
  boil::oout<<boil::endl;
  for(auto c : coefficients)
    boil::oout<<c<<" ";
  boil::oout<<boil::endl;

  std::vector<StencilPoint> stencil;
  int ii(0);
  for(auto s : stenpos) {
    stencil.push_back(
      StencilPoint(ii,evaluate_polynomial(1,coefficients,s),s));
    ii++;
  }

  /* first-order */
  real derivative_1st = coefficients[1];
  
  real difference_1st = first_order_first(stencil);

  flag = flag & fabs(derivative_1st-difference_1st)<boil::nano;
  boil::oout<<derivative_1st<<" "<<difference_1st<<" "<<flag<<boil::endl;

  /* second-order */
  for(auto & s : stencil) {
    s.val = evaluate_polynomial(2,coefficients,s.pos);
  } 
  real derivative_2nd = coefficients[1];
  real difference_2nd = second_order_first(stencil);

  flag = flag & fabs(derivative_2nd-difference_2nd)<boil::nano;
  boil::oout<<derivative_2nd<<" "<<difference_2nd<<" "<<flag<<boil::endl;

  /* third-order */
  for(auto & s : stencil) {
    s.val = evaluate_polynomial(3,coefficients,s.pos);
  } 
  real derivative_3rd = coefficients[1];
  real difference_3rd = third_order_first(stencil);

  flag = flag & fabs(derivative_3rd-difference_3rd)<boil::nano;
  boil::oout<<derivative_3rd<<" "<<difference_3rd<<" "<<flag<<boil::endl;

  /* fourth-order */
  for(auto & s : stencil) {
    s.val = evaluate_polynomial(4,coefficients,s.pos);
  } 
  real derivative_4th = coefficients[1];
  real difference_4th = fourth_order_first(stencil);

  flag = flag & fabs(derivative_4th-difference_4th)<boil::nano;
  boil::oout<<derivative_4th<<" "<<difference_4th<<" "<<flag<<boil::endl;


  return flag;

}

/******************************************************************************/
real Topology::evaluate_polynomial(const int order,
                                   const std::vector<real> & coefficients,
                                   const real x) {
/******************************************************************************/
  real val(0.);
  for(int i(0); i<=order; ++i) {
    val += coefficients[i]*pow(x,real(i)); 
  }

  return val;
}

/******************************************************************************/
real Topology::evaluate_polynomial_derivative(const int order,
                                       const std::vector<real> & coefficients,
                                       const real x) {
/******************************************************************************/
  real val(0.);
  for(int i(1); i<=order; ++i) {
    val += real(i)*coefficients[i]*pow(x,real(i-1)); 
  }

  return val;
}
