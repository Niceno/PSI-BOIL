/**
 ** calculates the disjoining pressure (p_l - p_v) 
 ** under the assumption that solid-liquid boundary is a plane
 ** with surface normal (0,0,1) and located at z = 0.
 **/
void calculate_disjoint_pressure(Scalar & dp, const Vector & bndclr,
                                 const Scalar & curv, const real sigma, 
                                 const real hamaker, const real delta0) {

  for_vijk(dp,i,j,k) {
    Comp m = Comp::k();
    const real ds = bndclr.dSz(m,i,j,k);
    const real clrb = std::min(std::max(bndclr[m][i][j][k  ]/ds,0.0),1.0);
    const real clrt = std::min(std::max(bndclr[m][i][j][k+1]/ds,0.0),1.0);

    if((clrb>0.5)&&(clrt<0.5)) { /* liquid to vapour transition */
      real coef;
      if(fabs(clrt-clrb)>boil::atto)
        coef = (0.5-clrb)/(clrt-clrb)-0.5;
      else
        coef = 0.0;
      real delta = coef*dp.dzc(k)+dp.zc(k);
      delta = std::max(delta,delta0);
      /* +sigma because curv < 0 for bubbles */
      dp[i][j][k] = +sigma*curv[i][j][k] - hamaker/pow(delta,3.0);
    } else {
      dp[i][j][k] = +sigma*curv[i][j][k];
    }
  }
  dp.exchange();
}

void calculate_disjoint_pressure_x(Scalar & dp, const Vector & bndclr,
                                   const Scalar & curv, const real sigma, 
                                   const real hamaker, const real delta0) {

  for_vijk(dp,i,j,k) {
    Comp m = Comp::i();
    const real ds = bndclr.dSx(m,i,j,k);
    const real clrb = std::min(std::max(bndclr[m][i  ][j][k]/ds,0.0),1.0);
    const real clrt = std::min(std::max(bndclr[m][i+1][j][k]/ds,0.0),1.0);

    if((clrb>0.5)&&(clrt<0.5)) { /* liquid to vapour transition */
      real coef;
      if(fabs(clrt-clrb)>boil::atto)
        coef = (0.5-clrb)/(clrt-clrb)-0.5;
      else
        coef = 0.0;
      real delta = coef*dp.dxc(i)+dp.xc(i);
      delta = std::max(delta,delta0);
      /* +sigma because curv < 0 for bubbles */
      dp[i][j][k] = +sigma*curv[i][j][k] - hamaker/pow(delta,3.0);
#if 0
      boil::oout << "Disj: " << i << " " << clrb << " "<<clrt<<" "<< dp.xc(i) << " "<< coef<<" "<< delta << " "<<delta0<< boil::endl;
#endif
    } else {
      dp[i][j][k] = +sigma*curv[i][j][k];
    }
  }
  dp.exchange();
}
