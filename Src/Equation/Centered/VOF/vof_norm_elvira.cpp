#include "vof.h"

/* Comp::i() == YZ-plane and 1 = y, 2 = z, 3 = x.
   Comp::j() == XZ-plane and 1 = x, 2 = z, 3 = y. 
   Comp::k() == XY-plane and 1 = x, 2 = y, 3 = z. 
   The 3-direction is not-considered and is identically zero. */
void VOF::norm_elvira(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Elvira method: Pilliod and Puckett, LBNL-40744 (1997), 2D
*         Results: nx, ny, nz  +  nalpha
*******************************************************************************/

  assert(mcomp_for_elvira != Comp::undefined());

  /* 2D algorithm */
  for_ijk(i,j,k) {
    norm_elvira_kernel(nx[i][j][k],ny[i][j][k],nz[i][j][k], i,j,k, sca);
    if(dom->ibody().off(i,j,k)) {
      nalpha[i][j][k] = ((sca[i][j][k]>=phisurf)-(sca[i][j][k]<phisurf))*boil::unreal;
    }
  }

  /* missing boundary conditions */
  insert_bc_norm_cc(sca);
  /* nalpha exchange inside extract alpha */
  extract_alpha_near_bnd(sca);

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-mx-my-mz", time->current_step());
  //exit(0);

  return;
}

void VOF::norm_elvira_kernel(real & nx_val, real & ny_val, real & nz_val,
                             const int i, const int j, const int k,
                             const Scalar & sca) {
  real valcc,valmc,valpc,valcm,valcp,valmm,valpm,valmp,valpp;

  if       (mcomp_for_elvira==Comp::i()) {
    valcc = sca[i][j  ][k  ];
    valmc = sca[i][j-1][k  ];
    valpc = sca[i][j+1][k  ];
    valcm = sca[i][j  ][k-1];
    valcp = sca[i][j  ][k+1];
    valmm = sca[i][j-1][k-1];
    valpm = sca[i][j+1][k-1];
    valmp = sca[i][j-1][k+1];
    valpp = sca[i][j+1][k+1];
  } else if(mcomp_for_elvira==Comp::j()) {
    valcc = sca[i  ][j][k  ];
    valmc = sca[i-1][j][k  ];
    valpc = sca[i+1][j][k  ];
    valcm = sca[i  ][j][k-1];
    valcp = sca[i  ][j][k+1];
    valmm = sca[i-1][j][k-1];
    valpm = sca[i+1][j][k-1];
    valmp = sca[i-1][j][k+1];
    valpp = sca[i+1][j][k+1];
  } else if(mcomp_for_elvira==Comp::k()) {
    valcc = sca[i  ][j  ][k];
    valmc = sca[i-1][j  ][k];
    valpc = sca[i+1][j  ][k];
    valcm = sca[i  ][j-1][k];
    valcp = sca[i  ][j+1][k];
    valmm = sca[i-1][j-1][k];
    valpm = sca[i+1][j-1][k];
    valmp = sca[i-1][j+1][k];
    valpp = sca[i+1][j+1][k];
  } else {
    boil::aout<<"Elvira direction not properly set! Exiting."<<boil::endl;
    exit(0);
  }

#if 1
  /* used only for cells in the vicinity of the interface
     to reduce computational cost */
  if(  (valcc-phisurf)*(valmc-phisurf)<=0.0
     ||(valcc-phisurf)*(valpc-phisurf)<=0.0
     ||(valcc-phisurf)*(valcm-phisurf)<=0.0
     ||(valcc-phisurf)*(valcp-phisurf)<=0.0
     ||(valcc-phisurf)*(valmm-phisurf)<=0.0
     ||(valcc-phisurf)*(valpm-phisurf)<=0.0
     ||(valcc-phisurf)*(valmp-phisurf)<=0.0
     ||(valcc-phisurf)*(valpp-phisurf)<=0.0) {
    norm_elvira_kernel_full(nx_val, ny_val, nz_val, mcomp_for_elvira,
                            i,j,k,valcc,valmc,valpc,valcm,valcp,
                                        valmm,valpm,valmp,valpp);
  /* otherwise, standard method is used */
  } else {
    norm_mixed_kernel(nx[i][j][k], ny[i][j][k], nz[i][j][k], i,j,k, sca);
    real scpval = sca[i][j][k];
    if(scpval==0.5) scpval += boil::pico;
    nalpha[i][j][k] = alpha_val(scpval,
                                nx[i][j][k],ny[i][j][k],nz[i][j][k]);
  }
#else
  norm_elvira_kernel_full(nx_val, ny_val, nz_val, mcomp_for_elvira,
                          i,j,k,valcc,valmc,valpc,valcm,valcp,
                                      valmm,valpm,valmp,valpp);
#endif

  return;
}

/* Comp::i() == YZ-plane and 1 = y, 2 = z, 3 = x.
   Comp::j() == XZ-plane and 1 = x, 2 = z, 3 = y. 
   Comp::k() == XY-plane and 1 = x, 2 = y, 3 = z. 
   The 3-direction is not-considered and is identically zero. */
void VOF::norm_elvira_kernel_full(real & nx_val, real & ny_val, real & nz_val,
                                  const Comp & m, const int i, const int j, const int k,
                                  real valcc,real valmc,real valpc,real valcm,
                                  real valcp,real valmm,real valpm,real valmp,real valpp) {

  /* step one: slope candidates */
  real m1U, m1D, m1C;
  real m2U, m2D, m2C;

  m1U = (valcp+valcc+valcm)-(valmp+valmc+valmm);
  m1D = (valpp+valpc+valpm)-(valcp+valcc+valcm);
  m1C = 0.5*((valpp+valpc+valpm)-(valmp+valmc+valmm));

  m2U = (valpc+valcc+valmc)-(valpm+valcm+valmm);
  m2D = (valpp+valcp+valmp)-(valpc+valcc+valmc);
  m2C = 0.5*((valpp+valcp+valmp)-(valpm+valcm+valmm));

  /* step one and a half: maximal difference in the other direction
     determines sign of the ordinate */
  real sgn1, sgn2;

  if(fabs(m1U)>fabs(m1D)) {
    if(fabs(m1U)>fabs(m1C)) {
      sgn1 = (m1U>=0.0)-(m1U<0.0);
    } else {
      sgn1 = (m1C>=0.0)-(m1C<0.0);
    }
  } else {
    if(fabs(m1D)>fabs(m1C)) {
      sgn1 = (m1D>=0.0)-(m1D<0.0);
    } else {
      sgn1 = (m1C>=0.0)-(m1C<0.0);
    }
  }

  if(fabs(m2U)>fabs(m2D)) {
    if(fabs(m2U)>fabs(m2C)) {
      sgn2 = (m2U>=0.0)-(m2U<0.0);
    } else {
      sgn2 = (m2C>=0.0)-(m2C<0.0);
    }
  } else {
    if(fabs(m2D)>fabs(m2C)) {
      sgn2 = (m2D>=0.0)-(m2D<0.0);
    } else {
      sgn2 = (m2C>=0.0)-(m2C<0.0);
    }
  }

  /* step two: normal vector candidates */
  real n11U, n21U, n31U;
  real n11D, n21D, n31D;
  real n11C, n21C, n31C;
  real n12U, n22U, n32U;
  real n12D, n22D, n32D;
  real n12C, n22C, n32C;

  normalize_elvira(m1U, sgn2, n11U, n21U, n31U);
  normalize_elvira(m1D, sgn2, n11D, n21D, n31D);
  normalize_elvira(m1C, sgn2, n11C, n21C, n31C);

  normalize_elvira(m2U, sgn1, n22U, n12U, n32U);
  normalize_elvira(m2D, sgn1, n22D, n12D, n32D);
  normalize_elvira(m2C, sgn1, n22C, n12C, n32C);

  /* step three: alpha candidates
   * note that alpha is valid for positive normal vector 
   * in the reduced cell (0.0,1.0)x(0.0,1.0)x(0.0,1.0) */
  real alp1U, alp1D, alp1C;
  real alp2U, alp2D, alp2C;

  /* for empty cells: */
  if       (valcc<boil::micro) {
    alp1U = alp1D = alp1C = alp2U = alp2D = alp2C = 0.0;
  /* for full cells: */
  } else if(valcc-1.>-boil::micro) {
    alp1U = fabs(n11U)+fabs(n21U)+fabs(n31U);
    alp1D = fabs(n11D)+fabs(n21D)+fabs(n31D);
    alp1C = fabs(n11C)+fabs(n21C)+fabs(n31C);
    alp2U = fabs(n12U)+fabs(n22U)+fabs(n32U); 
    alp2D = fabs(n12D)+fabs(n22D)+fabs(n32D); 
    alp2C = fabs(n12C)+fabs(n22C)+fabs(n32C);
  /* other cases */
  } else {
    alp1U = alpha_val(valcc,n11U,n21U,n31U);
    alp1D = alpha_val(valcc,n11D,n21D,n31D);
    alp1C = alpha_val(valcc,n11C,n21C,n31C);
    alp2U = alpha_val(valcc,n12U,n22U,n32U); 
    alp2D = alpha_val(valcc,n12D,n22D,n32D);  
    alp2C = alpha_val(valcc,n12C,n22C,n32C);
  }

  /* step four: for the 3x3 grid, evaluate the L2 error */
  real err1U, err1D, err1C, err2U, err2D, err2C;

  err1U = elvira_l2(alp1U,n11U,n21U,n31U,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  err1D = elvira_l2(alp1D,n11D,n21D,n31D,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  err1C = elvira_l2(alp1C,n11C,n21C,n31C,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);

  err2U = elvira_l2(alp2U,n12U,n22U,n32U,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  err2D = elvira_l2(alp2D,n12D,n22D,n32D,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  err2C = elvira_l2(alp2C,n12C,n22C,n32C,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);

  /* step five: choose the one, which minimises the error */
  std::vector<real> errs = { err1U, err1D, err1C, err2U, err2D, err2C };
  std::vector<real> n1s = { n11U, n11D, n11C, n12U, n12D, n12C };
  std::vector<real> n2s = { n21U, n21D, n21C, n22U, n22D, n22C };
  std::vector<real> n3s = { n31U, n31D, n31C, n32U, n32D, n32C };
  std::vector<real> alps = { alp1U, alp1D, alp1C, alp2U, alp2D, alp2C };
  int min_err_idx = std::min_element(errs.begin(),errs.end()) - errs.begin();

  if       (m==Comp::i()) {
    ny_val = n1s[min_err_idx];
    nz_val = n2s[min_err_idx];
    nx_val = n3s[min_err_idx];
  } else if(m==Comp::j()) {

#if 0
    if(i==3&&j==3&&k==34) {
    if(fabs(nx[i][j][k]-n1s[min_err_idx])>boil::pico&&adens[i][j][k]>boil::pico) {
      real theta = atan(phi.zc(k)/(phi.xc(i)+boil::atto));
      real nxtest = cos(theta);
      real nztest = sin(theta);
      boil::oout<<i<<" "<<j<<" "<<k<<" "<<phi[i][j][k]<<" "<<adens[i][j][k]<<" | "<<n1s[min_err_idx]<<" "<<n2s[min_err_idx]<<" "<<n3s[min_err_idx]<<" | "<<nx[i][j][k]<<" "<<nz[i][j][k]<<" "<<ny[i][j][k]<<" | "<<nxtest<<" "<<nztest<<boil::endl;
    }
#endif

    nx_val = n1s[min_err_idx];
    nz_val = n2s[min_err_idx];
    ny_val = n3s[min_err_idx];
  } else {
    nx_val = n1s[min_err_idx];
    ny_val = n2s[min_err_idx];
    nz_val = n3s[min_err_idx];
  }

  /* nalpha is also assigned */
#if 0
  real denom = fabs(n1s[min_err_idx])+fabs(n2s[min_err_idx])+fabs(n3s[min_err_idx]);
  nalpha[i][j][k] = denom*alps[min_err_idx]; 
  boil::oout<<i<<" "<<j<<" "<<k<<" | "<<alps[min_err_idx]<<" "<<alpha_val(valcc,nx_val,nz_val,ny_val)<<boil::endl;
#else
  nalpha[i][j][k] = alps[min_err_idx];
#endif

#if 0
  if((i==6)||(j==6)) {
    //boil::oout<<i<<" "<<j<<" | "<<nx[i][j][k]<<" "<<ny[i][j][k]<<" "<<nz[i][j][k]<<boil::endl;
    boil::oout<<i<<" "<<j<<" "<<valcc<<" "<<min_err_idx<<" | "<<nx[i][j][k]<<" "<<ny[i][j][k]<<" "<<nz[i][j][k]<<" | ";
    //if(i==10&&j==17)
      for(auto v : errs)
        boil::oout<<v<<" ";
    //else
      //for(auto v : nys)
        //boil::oout<<v<<" ";
    boil::oout<<boil::endl;
  }
#endif

  return;
}

void VOF::normalize_elvira(const real m, const real sig,
                           real & nn1, real & nn2, real & nn3) {

  nn1 = m;
  nn2 = sig;
  nn3 = 0.0;
  normalize(nn1,nn2,nn3);

  return;
}

#if 0
real VOF::elvira_l2(const int i, const int j, const int k,
                    const real alp, real & nnx, real & nny, real & nnz,
                    const real valcc,const real valmc,const real valpc,
                    const real valcm,const real valcp,const real valmm,
                    const real valpm,const real valmp,const real valpp) {
#else
real VOF::elvira_l2(const real alp,const real nn1,const real nn2,const real nn3,
                    const real valcc,const real valmc,const real valpc,
                    const real valcm,const real valcp,const real valmm,
                    const real valpm,const real valmp,const real valpp) {
#endif

  /* unnormalized alpha value */
  real alphaval = alp;

  if(!boil::realistic(alphaval)) {
    /* all extrapolated trivially */
    return (valmm-valcc)*(valmm-valcc)+(valmc-valcc)*(valmc-valcc) 
          +(valmp-valcc)*(valmp-valcc)+(valcm-valcc)*(valcm-valcc)
          +(valcp-valcc)*(valcp-valcc)+(valpm-valcc)*(valpm-valcc)
          +(valpc-valcc)*(valpc-valcc)+(valpp-valcc)*(valpp-valcc);
  }
          
  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nn1;
  real vn2 = -nn2;
  real vn3 = -nn3;

  /* normalized normal */
  real vm1 = fabs(vn1);
  real vm2 = fabs(vn2);
  real vm3 = fabs(vn3);

  real denom = vm1+vm2+vm3;

  /* step one: is the orientation correct? */
  real errmm_1 = valmm - ext_v(-0.5,-0.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 
  real errmc_1 = valmc - ext_v(-0.5, 0.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 
  real errmp_1 = valmp - ext_v(-0.5, 1.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 

  real errcm_1 = valcm - ext_v( 0.5,-0.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 
  /* errcc is by definition zero */
  real errcp_1 = valcp - ext_v( 0.5, 1.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 

  real errpm_1 = valpm - ext_v( 1.5,-0.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 
  real errpc_1 = valpc - ext_v( 1.5, 0.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 
  real errpp_1 = valpp - ext_v( 1.5, 1.5,0.5,vm1,vm2,vm3,vn1,vn2,vn3,denom,alphaval); 

  real toterr_1 = errmm_1*errmm_1+errmc_1*errmc_1
                + errmp_1*errmp_1+errcm_1*errcm_1
                + errcp_1*errcp_1+errpm_1*errpm_1
                + errpc_1*errpc_1+errpp_1*errpp_1;

#if 0 /* with the sgn determination, it should be always correct */
  /* step two: or is it inverted? */
  real errmm_2 = valmm - ext_v(-0.5,-0.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  real errmc_2 = valmc - ext_v(-0.5, 0.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  real errmp_2 = valmp - ext_v(-0.5, 1.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);

  real errcm_2 = valcm - ext_v( 0.5,-0.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  real errcp_2 = valcp - ext_v( 0.5, 1.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);

  real errpm_2 = valpm - ext_v( 1.5,-0.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  real errpc_2 = valpc - ext_v( 1.5, 0.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  real errpp_2 = valpp - ext_v( 1.5, 1.5,0.5,vm1,vm2,vm3,-vn1,-vn2,-vn3,denom,alphaval);
  
  real toterr_2 = errmm_2*errmm_2+errmc_2*errmc_2
                + errmp_2*errmp_2+errcm_2*errcm_2
                + errcp_2*errcp_2+errpm_2*errpm_2
                + errpc_2*errpc_2+errpp_2*errpp_2;
#endif

#if 0
  if((i==17&&j==9)||(j==17&&i==9)) {
  boil::oout<<toterr_1<<" "<<toterr_2<<" | "<<errmm_1<<" "<<errmc_1<<" "<<errmp_1<<" "<<errcm_1<<" "<<errcp_1<<" "<<errpm_1<<" "<<errpc_1<<" "<<errpp_1<<" | "
  //boil::oout<<toterr_2<<" | "<<valmm<<" "<<valmc<<" "<<valmp<<" "<<valcm<<" "<<valcp<<" "<<valpm<<" "<<valpc<<" "<<valpp<<" | "
            <<errmm_2<<" "<<errmc_2<<" "<<errmp_2<<" "<<errcm_2<<" "<<errcp_2<<" "<<errpm_2<<" "<<errpc_2<<" "<<errpp_2<<" | "<< boil::endl;
  }
#endif

#if 0
  real toterr;

  /* it was inverted after all */
  if(toterr_2<toterr_1) {
    toterr = toterr_2;
 
    /* don't forget to invert the normal vector */
    nnx = -nnx;
    nny = -nny;
    nnz = -nnz;
  } else {
    toterr = toterr_1;
  }

  return toterr;
#else
  return toterr_1;
#endif

}

/* the designations x,y,z are inconsequential, they merely refer to
   components 1,2,3 in the overall elvira call */
real VOF::ext_v(const real xp, const real yp, const real zp, 
                const real vv1, const real vv2, const real vv3,
                const real vn1, const real vn2, const real vn3, 
                const real denom, const real alp) {

  real xpos = xp;
  real ypos = yp;
  real zpos = zp;

  real vm1 = vv1;
  real vm2 = vv2;
  real vm3 = vv3;

  real alphaval = alp;

  /* mirror boundary cell to the normalized space */
  if(vn1<0) {
    xpos = 1.0-xpos;
  }
  if(vn2<0) {
    ypos = 1.0-ypos;
  }
  if(vn3<0) {
    zpos = 1.0-zpos;
  }

  /* now, we need to translate the coordinate system */
  /* so that x'(center of bnd cell) = (0.5,0.5,0.5)  */
  /* x' = x - x(center of bnd cell) + (0.5,0.5,0.5)  */
  /* this affects alpha value: */
  /* m dot x' = alpha + m dot [(0.5,0.5,0.5) - x(cbc)] = alpha' */
  alphaval += vm1*(0.5-xpos) + vm2*(0.5-ypos) + vm3*(0.5-zpos);

  /* normalized alpha value */
  alphaval /= denom; 
  vm1 /= denom;
  vm2 /= denom;
  vm3 /= denom;

  /* volume fraction */
  return calc_v(alphaval,vm1,vm2,vm3);

}
