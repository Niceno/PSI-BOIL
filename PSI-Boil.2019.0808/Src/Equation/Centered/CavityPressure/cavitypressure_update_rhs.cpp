#include "cavitypressure.h"

/***************************************************************************//**
*  Similar to Pressure but corrects for inactive gas cells 
*  and adds free surface effect. 
*******************************************************************************/
real CavityPressure::update_rhs() {

  /*-------------------------------------------------------+
  |  computation of sources, taking care of immersed body  |
  +-------------------------------------------------------*/
  real tot_dia = 0.0;

  for_ijk(i,j,k) {
    if(in_gas(i,j,k)) {
      //fnew[i][j][k] = 0.0;
      /* maybe helps with convergence? */
      fnew[i][j][k] = A.c[i][j][k]*Pcav(i,j,k);
    } else if(dom->ibody().off_p(i,j,k)) {
      fnew[i][j][k] = 0.0;
    } else {
      
      real a_w = dSx(Sign::neg(),i,j,k) * time->dti();
      real a_e = dSx(Sign::pos(),i,j,k) * time->dti();
      real a_s = dSy(Sign::neg(),i,j,k) * time->dti();
      real a_n = dSy(Sign::pos(),i,j,k) * time->dti();
      real a_b = dSz(Sign::neg(),i,j,k) * time->dti();
      real a_t = dSz(Sign::pos(),i,j,k) * time->dti();

      if(dom->ibody().nccells() > 0) {
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
    } /* liquid cell */
  } /* ijk */

  /*---------------------------+
  |  compute error and source  |
  +---------------------------*/
  real err = 0.0;
  real sou = 0.0;

  for_ijk(i,j,k) {
    if(!in_gas(i,j,k)) {
      err += fnew[i][j][k] * fnew[i][j][k];
      sou += fnew[i][j][k];
    }
  }

  boil::cart.sum_real( &err );
  boil::cart.sum_real( &sou );
  boil::cart.sum_real( &tot_dia );

  if(tot_dia > boil::pico) {
    err = sqrt(err) / tot_dia;
    boil::oout << "CavityPressure::update_rhs; err sou dia = "
               << err << " " << sou << " " << tot_dia << boil::endl;
  } else {
    err = 0.0;
    boil::oout << "CavityPressure::update_rhs; err sou dia = "
               << err << " " << sou << " " << tot_dia << boil::endl;
  }

  /*----------------------------------------------------------+
  |  add external sources, free surface and boundary effects  |
  +----------------------------------------------------------*/
  for_ijk(i,j,k) {
    fnew[i][j][k] += fext[i][j][k]
                   + fold[i][j][k]
                   + fbnd[i][j][k];
  }

  return err;

}
