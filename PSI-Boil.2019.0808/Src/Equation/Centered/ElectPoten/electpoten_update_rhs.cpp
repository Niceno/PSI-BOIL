#include "electpoten.h"

/***************************************************************************//**
*  Called just before solving the linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*  Local array "fnew" represents \f$ \{f\} \f$ in the above equation.
*
*  \f$ fnew = - \int_S {\sigma_e \bf u \otimes \bf B} \, dS\f
*
*******************************************************************************/
real ElectPoten::update_rhs() {

  update_vel();
  
  /*-------------------------------------------------------+
  |  computation of sources, taking care of immersed body  |
  +-------------------------------------------------------*/
  real tot_dia = 0.0;

  for_ijk(i,j,k) {

    real a_w = dSx(Sign::neg(),i,j,k); // dSx, dSy, dSz are always positive
    real a_e = dSx(Sign::pos(),i,j,k); // regardress neg() nor pos()
    real a_s = dSy(Sign::neg(),i,j,k); // dSx(neg)=dSx(pos)=dSx
    real a_n = dSy(Sign::pos(),i,j,k); 
    real a_b = dSz(Sign::neg(),i,j,k); 
    real a_t = dSz(Sign::pos(),i,j,k); 

    Comp mi = Comp::u();
    Comp mj = Comp::v();
    Comp mk = Comp::w();
    real bw, be, bs, bn, bb, bt;
    real ufc,vfc,wfc,bx,by,bz;

    /* west */
    ufc=uf[mi][i][j][k];
    vfc=vf[mi][i][j][k];
    wfc=wf[mi][i][j][k];
    valueFace(&B,mi,i,j,k,&bx,&by,&bz);
    bw = a_w * fluid()->sigma_e(mi,i,j,k)   * (vfc*bz - wfc*by);

    /* east */
    ufc=uf[mi][i+1][j][k];
    vfc=vf[mi][i+1][j][k];
    wfc=wf[mi][i+1][j][k];
    valueFace(&B,mi,i+1,j,k,&bx,&by,&bz);
    be = a_e * fluid()->sigma_e(mi,i+1,j,k) * (vfc*bz - wfc*by);

    /* south */
    ufc=uf[mj][i][j][k];
    vfc=vf[mj][i][j][k];
    wfc=wf[mj][i][j][k];
    valueFace(&B,mj,i,j,k,&bx,&by,&bz);
    bs = a_s * fluid()->sigma_e(mj,i,j,k)   * (wfc*bx - ufc*bz);

    /* north */
    ufc=uf[mj][i][j+1][k];
    vfc=vf[mj][i][j+1][k];
    wfc=wf[mj][i][j+1][k];
    valueFace(&B,mj,i,j+1,k,&bx,&by,&bz);
    bn = a_n * fluid()->sigma_e(mj,i,j+1,k) * (wfc*bx - ufc*bz);

    /* bottom */
    ufc=uf[mk][i][j][k];
    vfc=vf[mk][i][j][k];
    wfc=wf[mk][i][j][k];
    valueFace(&B,mk,i,j,k,&bx,&by,&bz);
    bb = a_b * fluid()->sigma_e(mk,i,j,k)   * (ufc*by - vfc*bx);

    /* top */
    ufc=uf[mk][i][j][k+1];
    vfc=vf[mk][i][j][k+1];
    wfc=wf[mk][i][j][k+1];
    valueFace(&B,mk,i,j,k+1,&bx,&by,&bz);
    bt = a_t * fluid()->sigma_e(mk,i,j,k+1) * (ufc*by - vfc*bx);

    fnew[i][j][k] = bw - be + bs - bn + bb - bt;

    tot_dia += sqrt(bw*bw + be*be + bs*bs + bn*bn + bb*bb + bt*bt);
#if 0
    if(i==5&&j==5)
    std::cout<<"update_rhs:dai= "<<k
             <<" "<<fnew[i][j][k]<<" "<<a_b<<" "
             <<uf[mk][i][j][k]<<" "<<uf[mk][i][j][k+1]<<"\n";
#endif
  }

  /*---------------------------+
  |  compute error and source  |
  +---------------------------*/
  real err = 0.0;
  real sou = 0.0;

  for_ijk(i,j,k) {
    err += fnew[i][j][k] * fnew[i][j][k];
    sou += fnew[i][j][k];                  
  }

  boil::cart.sum_real( &err );
  boil::cart.sum_real( &sou );
  boil::cart.sum_real( &tot_dia );

  if(tot_dia > boil::pico) {
    err = sqrt(err) / tot_dia;    
  } else {
    err = 0.0;
  }
  boil::oout << "@electpot_update_rhs; abs(RHS) sum(RHS) dia = " 
             << err << " " << sou << " " << tot_dia << boil::endl;

  //debug: if( time->current_step() % 100 == 0 ) 
  //  boil::plot->plot((*u),fnew, "uvw-fnew", time->current_step());
  

  /*--------------------------------------------+ 
  |  add external sources and boundary effects  |
  +--------------------------------------------*/
  for_ijk(i,j,k) {
    fnew[i][j][k] += fext[i][j][k] + fbnd[i][j][k];
  }

  return 0.0;
}
