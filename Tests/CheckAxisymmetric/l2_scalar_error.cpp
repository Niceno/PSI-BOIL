/******************************************************************************/
real l2_scalar_error(const Scalar & sca, const Scalar & scb) {
/***************************************************************************//**
 \brief Calculate the L2 difference between two scalar fields.
    output: L2
*******************************************************************************/

  real l2(-0.0);
  int cnt(0);
  for_vijk(sca,i,j,k) {
    //if(dom->ibody().off(i,j,k)) continue;
    //real sdiff = fabs(std::max(0.0,std::min(1.0,sca[i][j][k]))
    //                 -std::max(0.0,std::min(1.0,scb[i][j][k])));
    real sdiff = fabs(sca[i][j][k]-scb[i][j][k]);
    l2 += sdiff*sdiff;
    cnt++;
  }
  boil::cart.sum_real(&l2);
  boil::cart.sum_int(&cnt);

  return l2;
}
