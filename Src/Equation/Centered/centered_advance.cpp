#include "centered.h"

/******************************************************************************/
void Centered::advance() {

  OPR( conv_ts.N() );

  for_ijk(i,j,k)
    fnew[i][j][k] = fold[i][j][k]
                  + cnew[i][j][k] * conv_ts.N();
               // + fext[][][] ...             


  //debug: boil::plot->plot(fnew, "fnew", time->current_step());
  //debug: boil::plot->plot(cnew, "cnew", time->current_step());

  /* it has not been implemented for fluid/solid (and will probably never be) */
  assert( !solid() );


  for_ijk(i,j,k) {
    real dvi = 1.0 / dV(i,j,k);

    real ci = 1.0/fluid()->cp(i,j,k);
    real ri = 1.0/fluid()->rho(i,j,k);

    phi[i][j][k] = time->dt() * ri * ci * dvi * fnew[i][j][k];
  }

  /* refresh buffers */
  phi.exchange();
}

/*-----------------------------------------------------------------------------+
 '$Id: centered_advance.cpp,v 1.13 2009/07/24 11:21:31 niceno Exp $'/
+-----------------------------------------------------------------------------*/
