#include "enthalpyfd.h"
#include "def.h"

void EnthalpyFD::inertial(const Old old) {
  inertial(phi,old);
}
/***************************************************************************//**
*  Calculates inertial contribution to rhs
*******************************************************************************/
void EnthalpyFD::inertial(Scalar & sca, const Old old) {

  /*---------------------------+
  |  fold = rho * cp * T / dt  |
  +---------------------------*/
  real dti = time->dti();

  if(old==Old::no) {
    for_ijk(i,j,k) {
      if(dom->ibody().on(i,j,k)){
        /* phase change: xor indicates change of phase */
        if(cht.topo->above_interface_old(i,j,k) ^ cht.topo->above_interface(i,j,k)) {
          sca[i][j][k] = cht.Tint(i,j,k);     /* crude code */
        }
      }
    }
    sca.bnd_update();
    sca.exchange();
  }

  /* no transport in solid */
  if( !solid() ) 
    for_ijk(i,j,k) {
      real c,r;
      if(cht.above_interface(i,j,k,old)) {
        c = cht.cpl(i,j,k);
      } else {
        c = cht.cpv(i,j,k);
      }
      fold[i][j][k] = c * dV(i,j,k) * sca[i][j][k] * dti;
    }
  /* with transport in solid */
  else {
    for_ijk(i,j,k) {
      const real fV = dom->ibody().fV(i,j,k);
      const real cs = solid()->cp (i,j,k);
      real cf = cht.cpl(i,j,k);
      if(!cht.above_interface(i,j,k,old)) {
        cf = cht.cpv(i,j,k);
      }

      fold[i][j][k] = (cf*fV + cs*(1.0-fV)) * dV(i,j,k)
                    * sca[i][j][k] * dti;
    }
  }

  return;
}
