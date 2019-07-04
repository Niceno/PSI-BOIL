#include "vof.h"

void VOF::norm_elvira(Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Elvira method: Pilliod and Puckett, LBNL-40744 (1997), 2D x-y
*         Results: nx, ny, nz
*******************************************************************************/

  /* 2D algorithm, assumed x and y */
  for_ijk(i,j,k) {
    real valcc = sca[i  ][j  ][k];
    real valmc = sca[i-1][j  ][k];
    real valpc = sca[i+1][j  ][k];
    real valcm = sca[i  ][j-1][k];
    real valcp = sca[i  ][j+1][k];
    real valmm = sca[i-1][j-1][k];
    real valpm = sca[i+1][j-1][k];
    real valmp = sca[i-1][j+1][k];
    real valpp = sca[i+1][j+1][k];

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
      norm_elvira(i,j,k,valcc,valmc,valpc,valcm,valcp,
                        valmm,valpm,valmp,valpp);
    /* otherwise, centered-columns are used */
    } else {
      //nx[i][j][k] = ny[i][j][k] = nz[i][j][k] = 0.0;
      norm_cc(sca,i,j,k);
    }
  }

  /* missing boundary conditions */
  insert_bc_norm_cc(sca);

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-mx-my-mz", time->current_step());
  //exit(0);

  return;
}

void VOF::norm_elvira(int i,int j,int k,
                      real valcc,real valmc,real valpc,real valcm,
                      real valcp,real valmm,real valpm,real valmp,real valpp) {

  /* step one: slope candidates */
  real mxU, mxD, mxC;
  real myU, myD, myC;

  mxU = (valcp+valcc+valcm)-(valmp+valmc+valmm);
  mxD = (valpp+valpc+valpm)-(valcp+valcc+valcm);
  mxC = 0.5*((valpp+valpc+valpm)-(valmp+valmc+valmm));

  myU = (valpc+valcc+valmc)-(valpm+valcm+valmm);
  myD = (valpp+valcp+valmp)-(valpc+valcc+valmc);
  myC = 0.5*((valpp+valcp+valmp)-(valpm+valcm+valmm));

  /* step one and a half: maximal difference in the other direction
     determines sign of the ordinate */
  real sgnx, sgny;

  if(fabs(mxU)>fabs(mxD)) {
    if(fabs(mxU)>fabs(mxC)) {
      sgnx = (mxU>=0.0)-(mxU<0.0);
    } else {
      sgnx = (mxC>=0.0)-(mxC<0.0);
    }
  } else {
    if(fabs(mxD)>fabs(mxC)) {
      sgnx = (mxD>=0.0)-(mxD<0.0);
    } else {
      sgnx = (mxC>=0.0)-(mxC<0.0);
    }
  }

  if(fabs(myU)>fabs(myD)) {
    if(fabs(myU)>fabs(myC)) {
      sgny = (myU>=0.0)-(myU<0.0);
    } else {
      sgny = (myC>=0.0)-(myC<0.0);
    }
  } else {
    if(fabs(myD)>fabs(myC)) {
      sgny = (myD>=0.0)-(myD<0.0);
    } else {
      sgny = (myC>=0.0)-(myC<0.0);
    }
  }

  /* step two: normal vector candidates */
  real nxxU, nyxU, nzxU;
  real nxxD, nyxD, nzxD;
  real nxxC, nyxC, nzxC;
  real nxyU, nyyU, nzyU;
  real nxyD, nyyD, nzyD;
  real nxyC, nyyC, nzyC;

  normalize_elvira(mxU, sgny, nxxU, nyxU, nzxU);
  normalize_elvira(mxD, sgny, nxxD, nyxD, nzxD);
  normalize_elvira(mxC, sgny, nxxC, nyxC, nzxC);

  normalize_elvira(myU, sgnx, nyyU, nxyU, nzyU);
  normalize_elvira(myD, sgnx, nyyD, nxyD, nzyD);
  normalize_elvira(myC, sgnx, nyyC, nxyC, nzyC);

  /* step three: alpha candidates
   * note that alpha is valid for positive normal vector 
   * in the reduced cell (0.0,1.0)x(0.0,1.0)x(0.0,1.0) */
  real alpxU, alpxD, alpxC;
  real alpyU, alpyD, alpyC;

  /* for empty cells: */
  if       (valcc<boil::micro) {
    alpxU = alpxD = alpxC = alpyU = alpyD = alpyC = 0.0;
  /* for full cells: */
  } else if(valcc-1.>-boil::micro) {
    alpxU = fabs(nxxU)+fabs(nyxU)+fabs(nzxU);
    alpxD = fabs(nxxD)+fabs(nyxD)+fabs(nzxD);
    alpxC = fabs(nxxC)+fabs(nyxC)+fabs(nzxC);
    alpyU = fabs(nxyU)+fabs(nyyU)+fabs(nzyU); 
    alpyD = fabs(nxyD)+fabs(nyyD)+fabs(nzyD); 
    alpyC = fabs(nxyC)+fabs(nyyC)+fabs(nzyC);
  /* other cases */
  } else {
    alpxU = alpha_val(valcc,nxxU,nyxU,nzxU);
    alpxD = alpha_val(valcc,nxxD,nyxD,nzxD);
    alpxC = alpha_val(valcc,nxxC,nyxC,nzxC);
    alpyU = alpha_val(valcc,nxyU,nyyU,nzyU); 
    alpyD = alpha_val(valcc,nxyD,nyyD,nzyD);  
    alpyC = alpha_val(valcc,nxyC,nyyC,nzyC);
  }

  /* step four: for the 3x3 grid, evaluate the L2 error */
  real errxU, errxD, errxC, erryU, erryD, erryC;

  errxU = elvira_l2(alpxU,nxxU,nyxU,nzxU,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  errxD = elvira_l2(alpxD,nxxD,nyxD,nzxD,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  errxC = elvira_l2(alpxC,nxxC,nyxC,nzxC,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);

  erryU = elvira_l2(alpyU,nxyU,nyyU,nzyU,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  erryD = elvira_l2(alpyD,nxyD,nyyD,nzyD,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);
  erryC = elvira_l2(alpyC,nxyC,nyyC,nzyC,
                    valcc,valmc,valpc,valcm,
                    valcp,valmm,valpm,valmp,valpp);

  /* step five: choose the one, which minimises the error */
  std::vector<real> errs = { errxU, errxD, errxC, erryU, erryD, erryC };
  std::vector<real> nxs = { nxxU, nxxD, nxxC, nxyU, nxyD, nxyC };
  std::vector<real> nys = { nyxU, nyxD, nyxC, nyyU, nyyD, nyyC };
  std::vector<real> nzs = { nzxU, nzxD, nzxC, nzyU, nzyD, nzyC };
  //std::vector<real> alps = { alpxU, alpxD, alpxC, alpyU, alpyD, alpyC };
  int min_err_idx = std::min_element(errs.begin(),errs.end()) - errs.begin();

  nx[i][j][k] = nxs[min_err_idx];
  ny[i][j][k] = nys[min_err_idx];
  nz[i][j][k] = nzs[min_err_idx];

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
                           real & nnx, real & nny, real & nnz) {

  nnx = m;
  nny = sig;
  nnz = 0.0;
  normalize(nnx,nny,nnz);

  return;
}

#if 0
real VOF::elvira_l2(const int i, const int j, const int k,
                    const real alp, real & nnx, real & nny, real & nnz,
                    const real valcc,const real valmc,const real valpc,
                    const real valcm,const real valcp,const real valmm,
                    const real valpm,const real valmp,const real valpp) {
#else
real VOF::elvira_l2(const real alp,const real nnx,const real nny,const real nnz,
                    const real valcc,const real valmc,const real valpc,
                    const real valcm,const real valcp,const real valmm,
                    const real valpm,const real valmp,const real valpp) {
#endif

  /* unnormalized alpha value */
  real alphaval = alp;
          
  /* calculate vn1, vn2, vn3: normal vector at cell center */
  /* n points to the liquid */
  real vn1 = -nnx;
  real vn2 = -nny;
  real vn3 = -nnz;

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


void VOF::norm_cc(Scalar & sca,int i,int j,int k) {
    real nxX, nyX, nzX;
    nxX = copysign(1.0,+(sca[i+1][j][k]-sca[i-1][j][k]));
    nyX = 0.5 * ( (sca[i+1][j+1][k]+sca[i][j+1][k]+sca[i-1][j+1][k])
                - (sca[i+1][j-1][k]+sca[i][j-1][k]+sca[i-1][j-1][k]));
    nzX = 0.0;
    normalize(nxX,nyX,nzX);

    real nxY, nyY, nzY;
    nxY = 0.5 * ( (sca[i+1][j-1][k]+sca[i+1][j][k]+sca[i+1][j+1][k])
                - (sca[i-1][j-1][k]+sca[i-1][j][k]+sca[i-1][j+1][k]));
    nyY = copysign(1.0,+(sca[i][j+1][k]-sca[i][j-1][k]));
    nzY = 0.0;
    normalize(nxY,nyY,nzY);

    if(fabs(nxX)<fabs(nyY)) {
      nx[i][j][k]=nxY;
      ny[i][j][k]=nyY;
      nz[i][j][k]=nzY;
    } else {
      nx[i][j][k]=nxX;
      ny[i][j][k]=nyX;
      nz[i][j][k]=nzX;
    }

    return;
}
