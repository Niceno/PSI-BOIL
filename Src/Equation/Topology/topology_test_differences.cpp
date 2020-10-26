#include "topology.h"

#include <random>

/******************************************************************************/
bool Topology::test_differences(const int count) {
/******************************************************************************/
  std::mt19937 rng(123);
  std::uniform_real_distribution<real> distR1(1.,10.);
  std::uniform_real_distribution<real> distR2(-10.,10.);
  std::uniform_int_distribution<int> distI(0,1);

  bool flag(true);
  for(int i(0); i<count; ++i) {
 
    boil::oout<<i<<boil::endl;

    /* completely randomised stencil with random polynomial */
    flag = flag & test_differences({0.,(distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng),
                                       (distI(rng)>0 ? 1 : -1) * distR1(rng)},

                                   {distR2(rng),distR2(rng),distR2(rng),
                                    distR2(rng),distR2(rng)});
  }

  return flag;

}


/******************************************************************************/
bool Topology::test_differences(const std::vector<real> & stencil,
                                const std::vector<real> & coefficients) {
/***************************************************************************//**
*  \brief Test if differences work correctly using polynomials
*******************************************************************************/

  bool flag(true);

  for(auto s : stencil)
    boil::oout<<s<<" ";
  boil::oout<<boil::endl;
  for(auto c : coefficients)
    boil::oout<<c<<" ";
  boil::oout<<boil::endl;

  /* first-order */
  real derivative_1st = coefficients[1];
  real difference_1st = first_order_difference(stencil,
                        {evaluate_polynomial(1,coefficients,stencil[0]),
                         evaluate_polynomial(1,coefficients,stencil[1])});

  flag = flag & fabs(derivative_1st-difference_1st)<boil::nano;
  boil::oout<<derivative_1st<<" "<<difference_1st<<" "<<flag<<boil::endl;

  /* second-order */
  real derivative_2nd = coefficients[1];
  real difference_2nd = second_order_difference(stencil,
                        {evaluate_polynomial(2,coefficients,stencil[0]),
                         evaluate_polynomial(2,coefficients,stencil[1]),
                         evaluate_polynomial(2,coefficients,stencil[2])});

  flag = flag & fabs(derivative_2nd-difference_2nd)<boil::nano;
  boil::oout<<derivative_2nd<<" "<<difference_2nd<<" "<<flag<<boil::endl;

  /* third-order */
  real derivative_3rd = coefficients[1];
  real difference_3rd = third_order_difference(stencil,
                        {evaluate_polynomial(3,coefficients,stencil[0]),
                         evaluate_polynomial(3,coefficients,stencil[1]),
                         evaluate_polynomial(3,coefficients,stencil[2]),
                         evaluate_polynomial(3,coefficients,stencil[3])});

  flag = flag & fabs(derivative_3rd-difference_3rd)<boil::nano;
  boil::oout<<derivative_3rd<<" "<<difference_3rd<<" "<<flag<<boil::endl;

  /* fourth-order */
  real derivative_4th = coefficients[1];
  real difference_4th = fourth_order_difference(stencil,
                        {evaluate_polynomial(4,coefficients,stencil[0]),
                         evaluate_polynomial(4,coefficients,stencil[1]),
                         evaluate_polynomial(4,coefficients,stencil[2]),
                         evaluate_polynomial(4,coefficients,stencil[3]),
                         evaluate_polynomial(4,coefficients,stencil[4])});

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
