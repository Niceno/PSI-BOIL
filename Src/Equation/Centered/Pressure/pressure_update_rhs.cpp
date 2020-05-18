#include "pressure.h"

/***************************************************************************//**
*  Called just before solving the linear system of equations:
*  \f[
*      [A] \cdot \{ \phi \}^N = \{ f \}
*  \f] 
*  Local array "fnew" represents \f$ \{f\} \f$ in the above equation.
*
*  This implementation is quite different for Pressure than for other 
*  Centered variables. While in others it updates fnew with contribution 
*  from time-discretization of various terms (diffusion, convection, 
*  external), here it assembles only one contribution: divergence of the
*  velocity field.
*******************************************************************************/
real Pressure::update_rhs_pressure() {

  //u->exchange();
  
  /*-------------------------------------------------------+
  |  computation of sources, taking care of immersed body  |
  +-------------------------------------------------------*/
  real tot_dia = 0.0;

  for_ijk(i,j,k) {

    real a_w = dSx(Sign::neg(),i,j,k) * time->dti(); 
    real a_e = dSx(Sign::pos(),i,j,k) * time->dti(); 
    real a_s = dSy(Sign::neg(),i,j,k) * time->dti(); 
    real a_n = dSy(Sign::pos(),i,j,k) * time->dti(); 
    real a_b = dSz(Sign::neg(),i,j,k) * time->dti(); 
    real a_t = dSz(Sign::pos(),i,j,k) * time->dti(); 

    if(dom->ibody().nccells() > 0) 
      if( dom->ibody().on_p(i,j,k) ) {
        a_w *= dom->ibody().fSw(i,j,k);
        a_e *= dom->ibody().fSe(i,j,k);
        a_s *= dom->ibody().fSs(i,j,k);
        a_n *= dom->ibody().fSn(i,j,k);
        a_b *= dom->ibody().fSb(i,j,k);
        a_t *= dom->ibody().fSt(i,j,k);
      } 

    fnew[i][j][k] = 0.0;
    fnew[i][j][k] += a_w*(*u)[Comp::u()][i]  [j]  [k]  ; 
    fnew[i][j][k] -= a_e*(*u)[Comp::u()][i+1][j]  [k]  ;
    fnew[i][j][k] += a_s*(*u)[Comp::v()][i]  [j]  [k]  ; 
    fnew[i][j][k] -= a_n*(*u)[Comp::v()][i]  [j+1][k]  ;
    fnew[i][j][k] += a_b*(*u)[Comp::w()][i]  [j]  [k]  ; 
    fnew[i][j][k] -= a_t*(*u)[Comp::w()][i]  [j]  [k+1];

    if(dom->ibody().nccells() == 0) 
      tot_dia += ( fabs( a_w*(*u)[Comp::u()][i]  [j]  [k]  ) +
                   fabs( a_e*(*u)[Comp::u()][i+1][j]  [k]  ) +
                   fabs( a_s*(*u)[Comp::v()][i]  [j]  [k]  ) +
                   fabs( a_n*(*u)[Comp::v()][i]  [j+1][k]  ) +
                   fabs( a_b*(*u)[Comp::w()][i]  [j]  [k]  ) +
                   fabs( a_t*(*u)[Comp::w()][i]  [j]  [k+1]) );

    /*-------------------------------+
    |  a "touch" from immersed body  |
    +-------------------------------*/
    if(dom->ibody().nccells() > 0) {
      fnew[i][j][k] = 0.0;
      if( dom->ibody().on_p(i-1,j,k) ) {
        fnew[i][j][k] += a_w*(*u)[Comp::u()][i]  [j]  [k]  ; 
        tot_dia += fabs( a_w*(*u)[Comp::u()][i]  [j]  [k]  );
      }
      if( dom->ibody().on_p(i+1,j,k) ) {
        fnew[i][j][k] -= a_e*(*u)[Comp::u()][i+1][j]  [k]  ;
        tot_dia += fabs( a_e*(*u)[Comp::u()][i+1][j]  [k]  );
      }
      if( dom->ibody().on_p(i,j-1,k) ) {
        fnew[i][j][k] += a_s*(*u)[Comp::v()][i]  [j]  [k]  ; 
        tot_dia += fabs( a_s*(*u)[Comp::v()][i]  [j]  [k]  );
      }
      if( dom->ibody().on_p(i,j+1,k) ) {
        fnew[i][j][k] -= a_n*(*u)[Comp::v()][i]  [j+1][k]  ;
        tot_dia += fabs( a_n*(*u)[Comp::v()][i]  [j+1][k]  );
      }
      if( dom->ibody().on_p(i,j,k-1) ) {
        fnew[i][j][k] += a_b*(*u)[Comp::w()][i]  [j]  [k]  ; 
        tot_dia += fabs( a_b*(*u)[Comp::w()][i]  [j]  [k]  );
      }
      if( dom->ibody().on_p(i,j,k+1) ) {
        fnew[i][j][k] -= a_t*(*u)[Comp::w()][i]  [j]  [k+1];
        tot_dia += fabs( a_t*(*u)[Comp::w()][i]  [j]  [k+1]);
      }
    }
  }

  /*-------------------------------+
  |  a "touch" from immersed body  |
  +-------------------------------*/
  if(dom->ibody().nccells() > 0) 
    for_ijk(i,j,k) 
      if( dom->ibody().off_p(i,j,k) ) 
        fnew[i][j][k] = 0.0;

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
    boil::oout << "@pressure_update_rhs; err sou dia = " 
               << err << " " << sou << " " << tot_dia << boil::endl;
  } else {
    err = 0.0;
    boil::oout << "@pressure_update_rhs; err sou dia = " 
               << err << " " << sou << " " << tot_dia << boil::endl;
  }

  //debug: if( time->current_step() % 100 == 0 ) 
  //debug: boil::plot->plot(fnew, "p-fnew", time->current_step());

  /*-----------------------+ 
  |  add external sources  |
  +-----------------------*/
  for_ijk(i,j,k) {
    fnew[i][j][k] += fext[i][j][k];
  }

  return err;
}
