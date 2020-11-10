#include "enthalpyfd.h"
#include "def.h"

/***************************************************************************//**
*  \brief Interface for calling convection for new time step \f$\{C\}^{N}\f$.
*******************************************************************************/
void EnthalpyFD::convection() {
  convection(&cnew);
}

/***************************************************************************//**
*  \brief Computes convection term using the last available velocities. 
*
*  \param conv - array into which convection term will be stored,
*
*  If called from new_time_step(), it will create \f$ \{C\}^{N-1} \f$; 
*  if called from inner iteration loop (from SIMPLE algorithm) it will compute
*  \f$ \{C\}^{N} \f$. 
*******************************************************************************/
void EnthalpyFD::convection(Scalar * conv) {

  phi.exchange();

  *conv = 0.0;

#ifdef CNEW
  Old old = Old::no;
#else
  Old old = Old::yes;
#endif

  /*-------------------+
  |  calculate fluxes  |
  +-------------------*/
  for_m(m) {
    int ofx(0), ofy(0), ofz(0);
    if(m==Comp::u()) {
      if(!bflag_struct.ifull) continue;
      ofx++;
    }
    if(m==Comp::v()) {
      if(!bflag_struct.jfull) continue;
      ofy++;
    }
    if(m==Comp::w()) {
      if(!bflag_struct.kfull) continue;
      ofz++;
    }

    for_wvmijk(flux_liq,m,i,j,k) {
      /* velocity */
      real vel_liq = (*uliq)[m][i][j][k];
      real vel_gas = (*ugas)[m][i][j][k];

      /* value */
      real tval_liq = face_value(Sign::pos(),m,vel_liq,i,j,k,ofx,ofy,ofz,old);
      real tval_gas = face_value(Sign::neg(),m,vel_gas,i,j,k,ofx,ofy,ofz,old);

      /* flux */
      flux_liq[m][i][j][k] = vel_liq*tval_liq;
      flux_gas[m][i][j][k] = vel_gas*tval_gas;
    }
  }

  /*-------------------------+
  |  convection scalar (FD)  |
  +-------------------------*/
  /* conv = (-div(f) + T*div(u))*dV */
  for_ijk(i,j,k) {
    if(cht.above_interface(i,j,k,old)) {
      (*conv)[i][j][k] = flux_liq.divergence(i,j,k)
                       + phi[i][j][k]*uliq->divergence(i,j,k);
    } else {
      (*conv)[i][j][k] = flux_gas.divergence(i,j,k)
                       + phi[i][j][k]*ugas->divergence(i,j,k);
    }
    (*conv)[i][j][k] *= -dV(i,j,k);
  }

  /*----------------+
  |  immersed body  |
  +----------------*/
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) { 
        (*conv)[i][j][k] = 0.0;
      }
    }
  }


  /*----------------+
  |  heat capacity  |
  +----------------*/
  for_ijk(i,j,k) {
    real c;
    if(cht.above_interface(i,j,k,old)) {
      c = cht.cpl(i,j,k);
    } else {
      c = cht.cpv(i,j,k);
    }
    (*conv)[i][j][k] = c * (*conv)[i][j][k];
  }

  return;
} 
